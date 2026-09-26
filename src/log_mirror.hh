#pragma once

// Everything a program writes to standard output, mirrored to a log file in
// the working directory, as DIALS writes dials.find_spots.log and the rest.
// Standard error goes to the log too, in order: a log without the warning or
// the error that ended a run would mislead whoever reads it later.
//
// stdout and stderr are replaced with streams that write SYNCHRONOUSLY to the
// terminal and to the log. The usual built-in tee -- a pipe read by a thread --
// was not used, because a crash or a sanitizer report, written straight to the
// terminal just before the program aborts, can die in the pipe before the
// thread copies it out, and those are the times the output matters most. Here
// nothing is held anywhere but stdio's own line buffer, exactly as before.
//
// Header-only, so the spot finder, which does not link the pipeline's
// library, can use it too.

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <string>

#include <fcntl.h>
#include <unistd.h>

namespace mxi {

namespace log_detail {

struct Sink {
  int terminal = -1;
  int log = -1;
};

inline void write_all(int fd, const char *p, std::size_t n) {
  while (n > 0) {
    const ssize_t w = ::write(fd, p, n);
    if (w < 0) {
      if (errno == EINTR)
        continue;
      return;
    }
    p += w;
    n -= static_cast<std::size_t>(w);
  }
}

inline ssize_t sink_write(void *cookie, const char *buf, std::size_t n) {
  const Sink *s = static_cast<const Sink *>(cookie);
  // The terminal first: it is what someone is watching.
  write_all(s->terminal, buf, n);
  if (s->log >= 0)
    write_all(s->log, buf, n);
  return static_cast<ssize_t>(n);
}

#if defined(__APPLE__)
inline int sink_write_bsd(void *cookie, const char *buf, int n) {
  return static_cast<int>(sink_write(cookie, buf, static_cast<std::size_t>(n)));
}
inline std::FILE *open_sink(Sink *s) {
  return funopen(s, nullptr, sink_write_bsd, nullptr, nullptr);
}
#else
inline std::FILE *open_sink(Sink *s) {
  cookie_io_functions_t io{};
  io.write = sink_write;
  return fopencookie(s, "w", io);
}
#endif

} // namespace log_detail

//: True for a run that only prints help or a version. Such a run writes no
//: log, since it would otherwise overwrite the log of the last real run.
inline bool only_asks_for_help(int argc, char **argv) {
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "-h" || a == "--help" || a == "--version")
      return true;
  }
  return false;
}

//: From here on, stdout and stderr are mirrored into `path`, which is
//: truncated first. If the log cannot be opened, the program says so and
//: carries on without one: a missing log is no reason to lose a run.
inline void mirror_to_log(const char *path) {
  std::fflush(stdout);
  std::fflush(stderr);
  const int log = ::open(path, O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0644);
  if (log < 0) {
    std::fprintf(stderr,
                 "cannot write the log %s: %s; carrying on without it\n", path,
                 std::strerror(errno));
    return;
  }
  static log_detail::Sink out, err;
  out.terminal = STDOUT_FILENO;
  out.log = log;
  err.terminal = STDERR_FILENO;
  err.log = log;
  std::FILE *o = log_detail::open_sink(&out);
  std::FILE *e = log_detail::open_sink(&err);
  if (o == nullptr || e == nullptr) {
    std::fprintf(stderr, "cannot mirror output to %s; carrying on without it\n",
                 path);
    ::close(log);
    return;
  }
  // Line-buffered and unbuffered, as a terminal's stdout and stderr are, so
  // the two land in the log in the order they were written.
  std::setvbuf(o, nullptr, _IOLBF, 0);
  std::setvbuf(e, nullptr, _IONBF, 0);
  stdout = o;
  stderr = e;
}

} // namespace mxi
