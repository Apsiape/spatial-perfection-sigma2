Verification bundle, Version 1.0.0, November 2026

Run commands from verification/ with Python 3 on Windows, Linux or macOS:
  python review/check_good_root.py
  python review/check_tree_clock.py
  python compiler/exact_audit.py
  python compiler/minimal_compile.py
  python compiler/minimal_structural.py

The numerical tree suite requires NumPy. To avoid unnecessary BLAS threads,
set OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1 and MKL_NUM_THREADS=1.
The three compiler replication entry points use only the standard library.
Windows memory telemetry is optional; other platforms record null. Paths are relative to each script. The three entry points accept --output-dir DIR (or -o DIR), defaulting to compiler/generated/. They refuse the shipped compiler directory and toolkit subtree. Run minimal_compile.py and minimal_structural.py with the same output directory to audit a fresh compile; otherwise the structural audit copies the shipped read-only archives into its output directory.
Fresh stdout logs are in logs/. Exact result JSON files are in compiler/.

Provenance and modifications
review/: the three supplied review programs, unchanged. The toy tableau figure is generated from tableau(3,3)[3].
compiler/: independent compiler-replication programs and recorded archives.
The minimal_compile.py and minimal_structural.py source paths were changed
from machine-local locations to the bundled toolkit/ dependencies. No
mathematical rule was changed. The scripts generate adapted producer/checker
copies, removing their original hardcoded scenario size. The source hash in
minimal_results.json identifies the bundled original build_game.py.
The original tree review is numerical. The exact replication is a separate
implementation; its checker and producer dependencies are identified below.

Scope and fresh results
- Good roots: every prime below 4000; minimum good fraction 0.8738.
  The chord-bound subtest covers only primes below 400.
  The trigonometric test is floating point, not an exact all-prime proof.
- Numerical tree: (m,q,k)=(3,3,1),(3,5,2),(3,11,1),(4,3,1).
  All structural residuals agree to floating-point precision. Four diagonal
  isometry defects are -P_m, with mass 3/(4^(m+1)-1). Decoder is a toy front.
- exact_audit.py: 10 models (depths 0..4, q=2,3), 2000 exact checks,
  rational-pair arithmetic in Q(zeta_3), q=2 phase -1. Exactly four original
  rows fail; each is the top-depth source isometry row. The automaton is an
  explicit nonhalting pulse-sensitive front, not the general universal CA.
- minimal_compile.py: four bridges and one modular isometry; zero extra
  transported labels. 7 input polynomials, 34 words, 887 projection variables,
  748 projection rows, 27991 questions per party, at most 256 answers.
  Exact row contractions for cutoffs m=0..3 show only P_source failing.
  L=2,d=16; row event loss=t_m/46464; game loss=t_m/36404361907584.
- minimal_structural.py: exact FMS probabilities in Q(sqrt(2)); 2472 old
  forbidden events and 40 FMS PVMs. The adapted prior independent serialized
  audit checks 77750 bindings, 11652 column-graph rules, 748 row expansions,
  333 scalar pair contexts, spectator masks and character orthogonality.

The minimal compiler producer reuses the earlier row compiler. The structural
verifier reuses the earlier independent serialized checker. It is not a second
independent full game producer. No dense full-game matrix is allocated.
Neither the minimal instance nor finite numerical sampling proves the general
compiler, the full universal computation reduction, compactness, or hardness.
Those assertions rest on the written proofs and the supplied written review.
No proof-assistant certification or dimension-independent robustness is claimed.

SHA256SUMS lists every bundled regular file except itself, relative to this
folder. Fresh generated/ and __pycache__/ directories are excluded. Re-running the entry points writes only to the selected output directory, leaving the bundled JSON/ZIP and their hashes unchanged. Regenerate SHA256SUMS only after intentionally updating shipped files.

Written replication review
The independent compiler review completed with verdict HOLDS WITH REPAIRS.
It proves the general four-bridge compiler and its dimension-independent
completeness estimate in the normal bounded single-shift scope. Its per-gadget
verdicts distinguish existing general tools from the new assembly and estimate.
No required gadget fails. It does not audit the arithmetic upper bound or Bell
attainment, and it does not export a universal-decoder transition table. The
manuscript includes its general defect-weight calculation and all listed repairs.

Revision checks
Run python run_revision_checks.py to repeat all five suites, simulate the
non-Windows telemetry branch with WinDLL forbidden, test output-directory
rejection, compare mathematical JSON and check shipped-artifact immutability.
Fresh logs and summary are preserved in logs/revision-checks/. This is a
branch simulation, not a native Linux or macOS execution.
Run python review/generate_toy.py --output-dir generated/toy to reproduce
the paper tableau without overwriting the shipped LaTeX table.
