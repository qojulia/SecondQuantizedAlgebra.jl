(* Cross-package comparison benchmarks: sneg (Mathematica).

   Scenario definitions, canonical keys, fairness contract and output format
   live in BENCHMARKS.md.

   Two-level encoding: a single fermionic mode f. For one site the fermion
   algebra (f f\[Dagger] + f\[Dagger] f = 1, f^2 = 0) is isomorphic to the
   sigma algebra, and bosons commute with fermions, so all Jaynes-Cummings
   results are identical to the spin encoding. This is also sneg's idiomatic
   usage (quantum impurity models).

   sneg is eager: nc[] canonicalizes on construction, so operations are timed
   directly. Indexed sums with a symbolic bound and symbolic expectation
   values are not expressible, so those scenarios are absent.

   Usage:  wolframscript -f sneg_bench.wls
*)

Needs["sneg`sneg`"];

capSeconds = 30.0;

snegbosonoperators[b];
snegfermionoperators[f];
snegrealconstants[wc, wa, g, w, J, U];

(* -- Known-answer equivalence checks -------------------------------------- *)
checks = {
    {"a a\[Dagger] = 1 + a\[Dagger]a",
        Expand[nc[b[AN], b[CR]]] === Expand[1 + nc[b[CR], b[AN]]]},
    {"a a a\[Dagger] = a\[Dagger]a a + 2a",
        Expand[nc[b[AN], b[AN], b[CR]]] ===
            Expand[nc[b[CR], b[AN], b[AN]] + 2 b[AN]]},
    {"\[Sigma]\[Minus]\[Sigma]\[Plus] + \[Sigma]\[Plus]\[Sigma]\[Minus] = 1",
        Expand[nc[f[AN], f[CR]] + nc[f[CR], f[AN]]] === 1}
};
Print["Known-answer equivalence checks:"];
Scan[Print["  ", If[#[[2]], "ok  ", "FAIL  "], #[[1]]] &, checks];
If[! AllTrue[checks, #[[2]] &],
    Print["equivalence checks failed; aborting"]; Exit[1]];

(* -- Scenario thunks ------------------------------------------------------- *)
jcBuild[] := Expand[
    wc nc[b[CR], b[AN]] + wa nc[f[CR], f[AN]] +
        g (nc[b[CR], f[AN]] + nc[b[AN], f[CR]])];
Hjc = jcBuild[];

nested[H_, op_, d_] := Nest[Expand[komutator[H, #]] &, op, d];

jcPow[p_] := Nest[Expand[nc[#, Hjc]] &, Hjc, p - 1];

aad = Expand[nc[b[AN], b[CR]]];
fockReorder[n_] := Nest[Expand[nc[#, aad]] &, aad, n - 1];

bhBuildM[M_] := Expand[
    Sum[w nc[b[CR, k], b[AN, k]], {k, M}] +
        J Sum[nc[b[CR, k], b[AN, k + 1]] + nc[b[CR, k + 1], b[AN, k]], {k, M - 1}] +
        U/2 Sum[nc[b[CR, k], b[CR, k], b[AN, k], b[AN, k]], {k, M}]];
Hbh = bhBuildM[8];

scenarios = {
    {"jc_build", jcBuild},
    {"jc_H2", Function[{}, Expand[nc[Hjc, Hjc]]]},
    {"jc_heisenberg", Function[{}, Expand[komutator[Hjc, b[AN]]]]},
    {"jc_nested_d1", Function[{}, nested[Hjc, f[AN], 1]]},
    {"jc_nested_d2", Function[{}, nested[Hjc, f[AN], 2]]},
    {"jc_nested_d3", Function[{}, nested[Hjc, f[AN], 3]]},
    {"jc_nested_d4", Function[{}, nested[Hjc, f[AN], 4]]},
    {"jc_nested_d5", Function[{}, nested[Hjc, f[AN], 5]]},
    {"jc_nested_d6", Function[{}, nested[Hjc, f[AN], 6]]},
    {"jc_nested_d7", Function[{}, nested[Hjc, f[AN], 7]]},
    {"jc_pow_n2", Function[{}, jcPow[2]]},
    {"jc_pow_n3", Function[{}, jcPow[3]]},
    {"jc_pow_n4", Function[{}, jcPow[4]]},
    {"jc_pow_n5", Function[{}, jcPow[5]]},
    {"jc_pow_n6", Function[{}, jcPow[6]]},
    {"jc_pow_n7", Function[{}, jcPow[7]]},
    {"fock_reorder_n2", Function[{}, fockReorder[2]]},
    {"fock_reorder_n4", Function[{}, fockReorder[4]]},
    {"fock_reorder_n6", Function[{}, fockReorder[6]]},
    {"fock_reorder_n8", Function[{}, fockReorder[8]]},
    {"fock_reorder_n10", Function[{}, fockReorder[10]]},
    {"fock_reorder_n12", Function[{}, fockReorder[12]]},
    {"fock_reorder_n14", Function[{}, fockReorder[14]]},
    {"bh_chain_M2", Function[{}, bhBuildM[2]]},
    {"bh_chain_M4", Function[{}, bhBuildM[4]]},
    {"bh_chain_M16", Function[{}, bhBuildM[16]]},
    {"bh_chain_M32", Function[{}, bhBuildM[32]]},
    {"bose_hubbard_build", Function[{}, bhBuildM[8]]},
    {"bose_hubbard_H2", Function[{}, Expand[nc[Hbh, Hbh]]]}
};

(* -- Timing ----------------------------------------------------------------- *)
(* Minimum per-call time under a fixed wall-clock budget, matching the Julia
   and Python harnesses: calibrate a loop count so one batch runs at least
   batchTarget seconds (amortizing timer resolution), then keep timing batches
   until `budget` seconds have elapsed (or maxBatches batches collected) and
   report the fastest per-call time. maxBatches is the early-stop sample cap:
   fast scenarios stop after ~maxBatches x batchTarget seconds instead of
   burning the full budget, matching the Julia/Python harnesses. Same statistic
   (minimum) and same per-scenario effort as the other packages, rather than
   RepeatedTiming's mean. *)
minTiming[thunk_] := Module[
    {loops = 1, batchTarget = 0.1, budget = 10.0, maxBatches = 50,
     best = Infinity, elapsed = 0., batches = 0, dt},
    While[First[AbsoluteTiming[Do[thunk[], {loops}]]] < batchTarget && loops < 10^8,
        loops *= 2
    ];
    While[elapsed < budget && batches < maxBatches,
        dt = First[AbsoluteTiming[Do[thunk[], {loops}]]];
        best = Min[best, dt / loops];
        elapsed += dt;
        batches += 1
    ];
    best
];

results = <||>;
capped = <||>;
Do[
    {key, thunk} = sc;
    thunk[];  (* warm-up *)
    trial = First[AbsoluteTiming[thunk[]]];
    If[trial > capSeconds,
        capped[key] = trial;
        Print["  ", key, ": capped (", NumberForm[trial, {4, 1}], " s/op)"],
        t = minTiming[thunk];
        results[key] = <|"time_ns" -> t 10.^9|>;
        Print["  ", key, ": ", NumberForm[t 10.^3, {6, 3}], " ms"]
    ],
    {sc, scenarios}
];

(* -- Output ------------------------------------------------------------------ *)
out = <|
    "meta" -> <|
        "package" -> "sneg",
        "version" -> ToString[$SnegVersion],
        "language" -> "Mathematica",
        "runtime_version" -> ToString[NumberForm[$VersionNumber, {3, 1}]],
        "date" -> DateString[Today, "ISODate"]
    |>,
    "results" -> results,
    "capped" -> capped
|>;
outdir = FileNameJoin[{DirectoryName[$InputFileName], "results"}];
If[! DirectoryQ[outdir], CreateDirectory[outdir]];
path = FileNameJoin[{outdir, "sneg.json"}];
Export[path, out, "JSON"];
Print["wrote ", path];
