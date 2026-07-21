# 2026-07-21: Fixing the hardcoded Boost include path

## The short version

**The algorithm did not change.** Not one line of the estimator, the mailman
multiplication, the jackknife, or the random vector generation was touched. FAME
computes exactly what it computed before, and the example outputs are unchanged.

What changed is how the project *builds*: where the compiler looks for a Boost
header, how CMake locates Boost, what the README tells you to install, and what git
tracks. Four files, one of them a one-line edit repeated three times.

If you only read one section, read this one. Everything below is the reasoning.

---

## 1. What broke

A user tried to compile FAME and got:

```
fatal error: /usr/include/boost/random.hpp: No such file or directory
```

They traced it themselves and reported the cause: `src/FAME.cpp` had the location of
`boost/random.hpp` written into the source as an absolute path.

```cpp
#include "/usr/include/boost/random.hpp"
```

Their suggested fix was to change it to `#include "boost/random.hpp"` and let the
compiler search for it. They also suggested adding Boost to the README, which never
listed it as a prerequisite.

Both suggestions were right. We applied a slightly wider version of each, for reasons
in section 3.

## 2. What we found digging in

Three things worth knowing before reading the changes.

**The hardcoded path appears three times, not once.** Besides `src/FAME.cpp:22`, the
same line sits in `src/ge_and_jack.cpp:26` and `src/ge_and_jack_cov.cpp:23`. Those two
files are not listed in `CMakeLists.txt` and nothing in the repository references
them, so they never triggered the reported error. They are legacy code carrying the
same latent bug.

**Boost is barely used.** The entire dependency comes down to three lines in each
file, at `src/FAME.cpp:2316-2319`:

```cpp
boost::mt19937 seedr;
seedr.seed(std::time(0));
boost::normal_distribution<> dist(0,1);
boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > z_vec(seedr, dist);
```

These generate the standard normal random vectors that drive the randomized trace
estimator. All three types are header-only, which matters: FAME needs Boost's headers
at compile time but links against no Boost library. That is why the fix required no
linker changes.

**CMake never knew about Boost.** `CMakeLists.txt` only had
`include_directories(include/)`, which points at the bundled Eigen. There was no
`find_package(Boost ...)` anywhere. The build worked on the original author's machine
purely because the absolute path happened to be correct there.

## 3. Why the reported one-line fix was not enough

Changing the include to a relative form is necessary. It is not sufficient.

`#include <boost/random.hpp>` tells the compiler "search your include path for
this." If Boost is not *on* that path, the build fails again with a nearly identical
error. The compiler's default path covers `/usr/include`, so the fix works on a Linux
box where Boost came from the system package manager. It does not cover:

- macOS with Homebrew, where Boost lands in `/opt/homebrew/include`
- conda environments, where it lands under the environment prefix
- HPC clusters using environment modules, which is how most users of a genetics tool
  like this one will encounter it

So the one-line fix would have moved the failure from one group of users to another.
Teaching CMake to *find* Boost, and then hand its location to the compiler, is what
makes the fix hold everywhere. That is change 2 below.

## 4. What we changed

### Change 1: the include itself (3 files)

`src/FAME.cpp:22`, `src/ge_and_jack.cpp:26`, `src/ge_and_jack_cov.cpp:23`:

```diff
- #include "/usr/include/boost/random.hpp"
+ #include <boost/random.hpp>
```

Angle brackets rather than the reporter's quotes. Both work, and the practical
difference is small: quotes search relative to the including file first, then fall
back to the include path. Angle brackets are the conventional form for an external
library, and they match how the other external dependency is already included a few
lines above (`#include <Eigen/Dense>`). Consistency with the surrounding code won.

The two legacy files were fixed even though nothing compiles them. It is the same
one-line change, and leaving a known-broken path in the tree only sets a trap for
whoever wires those files into the build later. They were not otherwise modified and
not deleted; that is a maintainer's call, noted in section 6.

### Change 2: teach CMake to find Boost

`CMakeLists.txt`:

```diff
- cmake_minimum_required(VERSION 2.8.12)
  project(EM_C++)
- include_directories(include/)
+ cmake_minimum_required(VERSION 3.5)
+
+ # Boost is used only for its header-only random number generators, so no
+ # Boost libraries need to be linked.
+ find_package(Boost REQUIRED)
+ include_directories(include/ ${Boost_INCLUDE_DIRS})
```

Three things happen here.

`find_package(Boost REQUIRED)` asks CMake to locate Boost. It is called with no
`COMPONENTS` argument, which is correct precisely because the usage is header-only.
Asking for components would make CMake hunt for compiled `.so` files that FAME does
not need and does not link.

`${Boost_INCLUDE_DIRS}` passes whatever CMake found to the compiler as an `-I` flag.
This is the part that makes the include in change 1 resolve on machines where Boost
is not in `/usr/include`.

`REQUIRED` changes the *shape* of the failure when Boost is genuinely missing.
Before, you got a confusing compile error deep in a source file. Now you get a clear
message at configure time, before compilation starts:

```
CMake Error at .../FindPackageHandleStandardArgs.cmake:137 (message):
  Could NOT find Boost (missing: Boost_INCLUDE_DIR)
```

That tells a user what to install. The old error told them a file was missing from a
directory they had no reason to think about.

The `cmake_minimum_required` bump from 2.8.12 to 3.5 is a separate problem found
along the way. CMake 4.0 removed compatibility with anything below 3.5 and refuses to
configure such a project at all. Anyone on a current CMake could not build FAME
regardless of their Boost setup. The bump is the minimum that clears it.

One cosmetic note: on CMake 3.30 and newer you will see a `CMP0167` deprecation
warning about the `FindBoost` module. It is harmless and the build proceeds.
Silencing it would mean raising the minimum version further than needed and cutting
off users on older systems.

### Change 3: the README

The Prerequisites section never mentioned Boost, which is the reporter's second
point. It now lists Boost alongside g++, cmake, and make, states that only headers
are needed, and gives install commands for apt, yum, Homebrew, and conda. It also
documents the escape hatch for non-standard locations:

```
cmake .. -DBOOST_ROOT=/path/to/boost
```

The stated cmake version was raised to match `CMakeLists.txt`.

While there, we added the Python packages the pipeline needs. `pipeline/*.py` import
numpy, pandas, scikit-learn, statsmodels, and bed_reader, and the README listed none
of them. Same class of problem as the Boost omission: a dependency you only discover
by hitting an error.

### Change 4: stop tracking the build directory

`build/` was committed to git: 32 files, 26 MB, including a compiled `GENIE_GxG`
binary and a `CMakeCache.txt` with paths baked in from the original author's cluster
account (`/u/scratch/b/boyang19/...`). A stale cache like that actively breaks a
rebuild, which is why the README opened its compile instructions by telling users to
delete the folder first.

We removed it from git tracking with `git rm -r --cached build/` (the files stay on
disk, they are just no longer version-controlled) and added a top-level `.gitignore`
containing `build/`.

With that done, the README line instructing users to delete the folder was removed.
It existed only to work around the committed artifact. A fresh clone no longer has a
`build/` directory to delete.

## 5. What we verified, and what we could not

Being precise here matters, because the verification was partial and you should know
which half.

**Confirmed on the development machine (macOS):**

1. No hardcoded paths remain anywhere in `src/` or `include/`.
2. `find_package(Boost REQUIRED)` locates Boost and passes it through. Checked by
   pointing the build at a stub Boost tree and reading the generated
   `flags.make`: both `include/` and the Boost directory appear as `-I` flags on the
   compile line. This is the mechanism the whole fix depends on.
3. `-DBOOST_ROOT=/path/to/boost` works exactly as the new README documents.
4. Missing Boost now fails at configure time with `Could NOT find Boost`, not with
   the reporter's compile-time error.
5. **The Boost include is no longer what blocks compilation.** With Boost reachable,
   the compiler now reads past line 22 and stops later, at line 8, on a different
   and pre-existing issue (see below). The reported bug is closed.

**Not verified:** the finished binary and the two bundled examples were not run.
The development machine is macOS with Apple clang, which cannot build FAME for
reasons that predate this change and are unrelated to Boost: `src/FAME.cpp:8`
includes `<bits/stdc++.h>`, a GCC-only header that clang does not ship, and
`CMAKE_CXX_FLAGS` passes `-lrt`, a library that does not exist on macOS. The README
scopes FAME to Linux, so these are expected. No Boost installation was present on
the machine either.

Because no numerical code changed, the example outputs cannot have shifted. To
confirm the full path on a Linux machine:

```bash
rm -rf build && mkdir build && cd build && cmake .. && make
cd ../example/single-bin && sh test.sh
cd ../multi-bin && sh test.sh
```

Compare against the committed reference values rather than diffing for exact
equality. FAME seeds its generator with `std::time(0)` (`src/FAME.cpp:2317`), so it
is already non-deterministic between runs by design. Results should land within
roughly one standard error of:

| Example | Expected |
|---------|----------|
| single-bin | `sigma^2_0` ~0.579 (se 0.181), `sigma^2_1` ~-0.008 (se 0.011), `sigma^2_2` ~0.406 (se 0.107) |
| multi-bin | `sigma^2_0` ~0.286 (se 0.116), `sigma^2_1` ~0.278 (se 0.101), `sigma^2_2` ~-0.007 (se 0.007), `sigma^2_3` ~0.420 (se 0.104) |

The `.results` files are tracked in git, so restore them with `git checkout` after
testing.

## 6. What we deliberately left alone

**Removing the Boost dependency entirely.** The three lines FAME uses map directly
onto C++11's `<random>` (`std::mt19937`, `std::normal_distribution`). The header is
even sitting commented out at `src/FAME.cpp:6`, and `-std=c++0x` is already set. That
would delete the prerequisite instead of documenting it. We did not do it because it
changes the exact sequence of random draws. The estimator is stochastic and
time-seeded, so results would stay statistically equivalent, but that is a
behavior-affecting change and belongs in its own decision, not smuggled into a build
fix.

**Deleting `src/ge_and_jack.cpp` and `src/ge_and_jack_cov.cpp`.** Both are absent
from `CMakeLists.txt` and unreferenced by every script in the repository. They look
like dead code and are probably worth removing, but that is a judgment only the
maintainers can make.

**`-lrt` and `<bits/stdc++.h>`.** Both are GCC and Linux specific, and both will
break a native macOS or clang build no matter what happens with Boost. The README
scopes the project to Linux, so they were left in place. Worth knowing that the
`brew install boost` line now in the README only helps macOS users who are already
building with a real GCC rather than Apple's clang.

## 7. Files touched

| File | Change |
|------|--------|
| `src/FAME.cpp` | Line 22: include path |
| `src/ge_and_jack.cpp` | Line 26: include path (file not built) |
| `src/ge_and_jack_cov.cpp` | Line 23: include path (file not built) |
| `CMakeLists.txt` | `find_package(Boost)`, include dirs, cmake minimum 2.8.12 to 3.5 |
| `README.md` | Boost and Python prerequisites, install commands, removed stale build-folder instruction |
| `.gitignore` | New file, ignores `build/` |
| `build/` | 32 files untracked (kept on disk) |
