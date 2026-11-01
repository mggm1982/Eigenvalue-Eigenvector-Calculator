EIGENVALUE/EIGENVECTOR CALCULATOR TEST DATA

Each file in inputs/ is a separate stdin test case.
Run one file at a time, e.g. in PowerShell:
  Get-Content .\inputs\3x3_01_distinct.txt | .\program.exe
Or in Command Prompt:
  program.exe < inputs\3x3_01_distinct.txt

Expected eigenvalues and sample eigenvectors are in expected_results.csv.
Eigenvectors can be multiplied by any nonzero scalar; different scaled answers may be correct.

WARNING: Some cases deliberately expose limitations of the supplied code:
- 2x2 eigenvector2x2 can divide by zero or choose an incorrect branch.
- Complex 2x2 eigenvalues cause uninitialized y1/y2 to be used.
- 3x3 cubicroot sign-change scanning can miss repeated eigenvalues.
- Numerical approximation may cause an incorrect rank or missing eigenvectors.
- The cubic search interval is fixed at -100 to 100.
