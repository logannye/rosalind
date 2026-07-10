/* tslint:disable */
/* eslint-disable */

/**
 * Incremental BLAKE3 hashing for large dropped artifacts. JavaScript should
 * feed fixed-size `File.slice()` chunks rather than materializing the whole file.
 */
export class Blake3Hasher {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Return the current digest without consuming the hasher.
     */
    finalize_hex(): string;
    /**
     * Create an empty hasher.
     */
    constructor();
    /**
     * Add one artifact chunk.
     */
    update(bytes: Uint8Array): void;
}

/**
 * Diff two receipts using the same cause/effect/noise localizer as the CLI.
 */
export function diff_receipts(a: string, b: string): string;

/**
 * Evaluate the shared trust model after JavaScript content-matches artifacts and
 * optionally supplies a reproduction certificate. `artifact_status` is one of
 * `not-checked`, `complete`, or `incomplete`.
 */
export function evaluate_trust(receipt_json: string, artifact_status: string, certificate_json: string): string;

/**
 * Inspect claim identity and resource/reproduction metadata for one receipt.
 */
export function inspect_receipt(json: string): string;

/**
 * Check receipt self-integrity. This intentionally does not claim that artifact
 * files were re-hashed or that the computation was reproduced.
 */
export function verify(json: string): string;

/**
 * Walk newline-delimited canonical receipts as a provenance DAG. Canonical
 * receipts contain no literal newlines, making JSONL a dependency-free bridge.
 */
export function walk_receipt_chain(receipts_jsonl: string): string;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly verify: (a: number, b: number) => [number, number];
    readonly inspect_receipt: (a: number, b: number) => [number, number];
    readonly evaluate_trust: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number];
    readonly diff_receipts: (a: number, b: number, c: number, d: number) => [number, number];
    readonly walk_receipt_chain: (a: number, b: number) => [number, number];
    readonly __wbg_blake3hasher_free: (a: number, b: number) => void;
    readonly blake3hasher_new: () => number;
    readonly blake3hasher_update: (a: number, b: number, c: number) => void;
    readonly blake3hasher_finalize_hex: (a: number) => [number, number];
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
