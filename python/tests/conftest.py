import sys
from pathlib import Path

# fixtures.py sits beside the tests and is imported by name rather than as part
# of the package: it generates test data and has no business shipping.
sys.path.insert(0, str(Path(__file__).parent))
