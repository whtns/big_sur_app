# BigSuR Performance + WASM Implementation Plan

Date: 2026-03-10
Owner: BigSuR engineering
Status: Proposed

## Goal
Improve perceived and actual performance for correlation workflows and large interactive visualizations, while using WebAssembly (WASM) only where it creates measurable wins.

## Decision Summary
- Prioritize low-risk backend and payload optimizations first.
- Run one narrow WASM pilot for client-side correlation network interaction.
- Keep Python/native path as default fallback until pilot proves ROI.

## Performance Targets
- Correlation tab open-to-first-render (warm cache): p95 <= 2.0s
- Correlation tab open-to-first-render (cold cache): p95 <= 8.0s
- Correlation network interaction (threshold slider, community switch): p95 <= 300ms
- Correlation task completion time on representative benchmark dataset: >= 25% improvement vs baseline OR equivalent UX gain via progressive rendering
- Browser main-thread long tasks (>50ms) during interaction: reduce by >= 40%

## Scope
In scope:
- Correlation computation/visualization path
- Dash callback payload and rendering efficiency
- Caching/read-write path improvements that impact interactive latency
- One Rust+WASM pilot for a single client-side hotspot

Out of scope:
- Full Scanpy/BigSur rewrite to WASM
- Replacing Celery architecture
- Rewriting all visualizations away from Plotly

## Phase 0: Benchmark and Instrumentation
### Ticket PERF-001: Add server-side timing instrumentation
Files:
- src/correlations/correlation_callbacks.py
- src/correlations/correlation_functions.py
- src/tasks/tasks.py
- src/processing/processing_functions.py

Tasks:
- Add structured timing logs around callback start/end and key compute blocks.
- Add payload size logging for dcc.Store objects and generated figure JSON.
- Add task-state timestamps for correlation jobs (queued/start/progress/success/failure).

Acceptance criteria:
- Every correlation request emits a request_id/session_id and stage timings.
- Logs include figure size (bytes) and node/edge counts for rendered networks.
- A one-command grep can extract p50/p95 by stage from logs.

### Ticket PERF-002: Add browser-side performance capture
Files:
- src/assets/jquery-3.4.1.min.js (do not edit)
- src/assets/MiCV.css (optional style hooks)
- new: src/assets/perf.js

Tasks:
- Add lightweight browser performance markers for tab load and interaction events.
- Track frame drops/long tasks around network redraw.
- Emit metrics to console and optional endpoint for local profiling runs.

Acceptance criteria:
- For each interaction, app records time-to-visual-update and long-task count.
- Instrumentation can be disabled with an env/config flag.

## Phase 1: No-WASM High-ROI Optimizations
### Ticket PERF-010: Remove dense matrix expansion in graph build path
Files:
- src/correlations/correlation_functions.py

Current issue:
- create_graph_from_mcPCCs calls todense then from_numpy_array, causing memory and CPU blowups for larger matrices.

Tasks:
- Build graph directly from sparse COO/CSR nonzero entries.
- Preserve edge weights and gene mapping exactly.
- Add guards for pathological sparse dimensions.

Acceptance criteria:
- No todense usage in correlation graph creation path.
- Graph node/edge counts match current behavior on test fixtures.
- Peak memory during graph creation reduced by >= 30% on benchmark.

### Ticket PERF-011: Reduce Plotly trace explosion for edges
Files:
- src/correlations/correlation_plotting.py

Current issue:
- One trace per edge creates heavy figure JSON and browser rendering overhead.

Tasks:
- Batch edges into minimal trace sets (for example grouped by sign/type/width bins).
- Preserve hover and legend behavior where needed.
- Keep visual parity for community/network views.

Acceptance criteria:
- Figure trace count reduced by >= 70% for representative large graph.
- Interaction latency improves by >= 30% on benchmark.
- Visual output validated against baseline snapshots.

### Ticket PERF-012: Shrink Dash store payloads and avoid repeated heavy serialization
Files:
- src/correlations/correlation_callbacks.py

Tasks:
- Store compact identifiers/paths instead of large dict payloads where possible.
- Avoid repeated DataFrame serialization/deserialization in hot callbacks.
- Cache immutable derived artifacts per session.

Acceptance criteria:
- dcc.Store payload size reduced by >= 50% in correlation flow.
- Callback CPU time for update_correlation_visuals reduced measurably (>= 20%).

### Ticket PERF-013: Improve polling strategy for background tasks
Files:
- src/correlations/correlation_components.py
- src/correlations/correlation_callbacks.py

Tasks:
- Replace fixed 500ms polling with adaptive backoff (fast initially, slower later).
- Stop polling immediately on terminal task states.

Acceptance criteria:
- Poll requests per task reduced by >= 40% with no UX regression.
- Time-to-status-update remains <= 1.0s p95 during active compute.

### Ticket PERF-014: Reduce unnecessary cache rewrites for AnnData
Files:
- src/helper_functions.py
- src/tasks/tasks.py

Tasks:
- Audit cache_adata write calls in hot paths and skip no-op rewrites.
- Separate metadata updates from full-object write where possible.
- Limit dense materialization background tasks to required layers only.

Acceptance criteria:
- Full zarr rewrite count per correlation interaction reduced by >= 50%.
- No data integrity regressions in session cache behavior.

### Ticket PERF-015: Restore joblib parallelism inside Celery forked workers
Files:
- src/tasks/tasks.py

Current issue:
- Celery uses a forked process pool (billiard). joblib's default Loky backend cannot spawn
  new processes inside a forked worker, so it silently falls back to n_jobs=1 for every
  call to calculate_correlations. Each correlation task (~2969 genes) runs fully single-threaded
  despite available cores, taking ~60s instead of <30s.

Tasks:
- Wrap calculate_correlations call with joblib.parallel_backend('threading', n_jobs=-1).
- Verify warning no longer appears in celery-worker.log.
- Benchmark task duration before/after on representative dataset.

Acceptance criteria:
- No "Loky-backed parallel loops cannot be called in a multiprocessing" warnings in logs.
- compute_correlations_task duration reduced by >= 25% vs baseline on benchmark dataset.
- Scientific output (mcPCCs, p-values) unchanged vs baseline reference vectors.

## Phase 2: WASM Pilot
### Ticket WASM-100: Define WASM pilot contract and fallback behavior
Files:
- new: docs/wasm_pilot_contract.md
- src/correlations/correlation_callbacks.py
- new: src/assets/wasm_loader.js

Tasks:
- Define exactly one operation to move client-side (recommended: edge threshold/filter + layout preprocessing for displayed subgraph).
- Define input/output schema and deterministic behavior.
- Add feature flag and fallback to existing Python/JS path.

Acceptance criteria:
- Feature can be toggled at runtime without restart.
- Fallback path produces equivalent output for same inputs.

### Ticket WASM-101: Implement Rust WASM module for selected hotspot
Files:
- new: wasm/correlation-kernel/Cargo.toml
- new: wasm/correlation-kernel/src/lib.rs
- new: wasm/correlation-kernel/README.md

Tasks:
- Implement chosen kernel (for example thresholding, edge bucketing, optional layout prep math).
- Use stable serialization format (typed arrays/compact JSON).
- Add deterministic tests in Rust for correctness.

Acceptance criteria:
- WASM kernel passes all correctness tests vs Python reference vectors.
- Per-operation runtime in browser is >= 2x faster than baseline JS/Python-client path on benchmark input sizes.

### Ticket WASM-102: Integrate WASM with Dash frontend callbacks
Files:
- src/correlations/correlation_components.py
- src/correlations/correlation_callbacks.py
- new: src/assets/correlation_wasm_bridge.js

Tasks:
- Wire WASM execution into correlation interaction flow.
- Keep UI responsive with non-blocking execution (worker or async chunking).
- Add user-visible fallback messaging for unsupported browsers/errors.

Acceptance criteria:
- No hard failures when WASM fails to load; auto-fallback works.
- Interaction p95 meets target (<= 300ms) on benchmark scenarios where WASM is enabled.

## Phase 3: Validation and Rollout
### Ticket PERF-200: Regression suite for performance-sensitive paths
Files:
- new: tests/perf/test_correlation_perf.py
- new: tests/perf/fixtures/README.md

Tasks:
- Add repeatable benchmark harness for representative datasets.
- Validate correctness parity for graphs/community outputs across optimized paths.

Acceptance criteria:
- CI/local benchmark command produces comparable report (baseline vs current).
- Correctness checks pass for all benchmark fixtures.

### Ticket PERF-201: Rollout checklist and observability
Files:
- README.md
- docs/troubleshooting.rst
- new: docs/performance_runbook.md

Tasks:
- Document feature flags, fallback paths, and profiling steps.
- Add runbook for diagnosing regressions and disabling WASM path quickly.

Acceptance criteria:
- Operators can enable/disable WASM without code changes.
- Runbook includes clear rollback steps and known failure modes.

## Milestones and Sequence
1. M1 Baseline ready: PERF-001, PERF-002
2. M2 No-WASM wins shipped: PERF-010..014
3. M3 WASM pilot integrated behind flag: WASM-100..102
4. M4 Performance regression harness + docs: PERF-200..201

## Risks and Mitigations
- Risk: WASM integration complexity exceeds value.
  - Mitigation: strict pilot scope, hard ROI gate before expansion.
- Risk: Optimizations alter scientific outputs.
  - Mitigation: parity tests against current output for nodes/edges/communities and ranking metrics.
- Risk: Browser variability (Safari/WebKit behavior, memory limits).
  - Mitigation: feature detection, fallback path, cross-browser smoke tests.

## ROI Gate for Continuing WASM
Proceed beyond pilot only if all are true:
- Interaction p95 improves by >= 30%
- Browser long tasks reduce by >= 40%
- No correctness regressions in parity suite
- Operational burden (build/deploy/debug) is acceptable to team

## Suggested Task Breakdown for Sprint Planning
- Sprint A: PERF-001, PERF-002, PERF-010
- Sprint B: PERF-011, PERF-012, PERF-013
- Sprint C: PERF-014, WASM-100, WASM-101 (partial)
- Sprint D: WASM-102, PERF-200, PERF-201

## Notes for Implementation
- Keep scientific compute server-side unless benchmark proves clear client-side gain.
- Prefer sparse-first data structures end-to-end in correlation workflow.
- Favor reducing bytes and trace count before adding new runtime complexity.
