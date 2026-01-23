/* tslint:disable */
/* eslint-disable */

export class MappingResult {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    bam_data: Uint8Array;
}

export class SingWebEngine {
    free(): void;
    [Symbol.dispose](): void;
    export_index(): Uint8Array;
    get_avg_coverage(): number;
    get_genome_size(): number;
    get_mapped_bases(): number;
    get_mapped_reads(): number;
    get_total_reads(): number;
    import_index(data: Uint8Array): void;
    map_reads_chunk(reads_chunk: Uint8Array): MappingResult;
    constructor();
    process_reference_chunk(chunk_data: Uint8Array, _chunk_id: number): void;
    run_mapping_chunk(ref_index_ptr: number, read1_chunk: Uint8Array, read2_chunk: Uint8Array | null | undefined, output_format: string, sort_output: boolean, write_header: boolean): MappingResult;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_mappingresult_free: (a: number, b: number) => void;
    readonly __wbg_get_mappingresult_bam_data: (a: number) => [number, number];
    readonly __wbg_set_mappingresult_bam_data: (a: number, b: number, c: number) => void;
    readonly __wbg_singwebengine_free: (a: number, b: number) => void;
    readonly singwebengine_new: () => number;
    readonly singwebengine_process_reference_chunk: (a: number, b: number, c: number, d: number) => [number, number];
    readonly singwebengine_export_index: (a: number) => [number, number, number, number];
    readonly singwebengine_import_index: (a: number, b: number, c: number) => [number, number];
    readonly singwebengine_run_mapping_chunk: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number, i: number, j: number) => [number, number, number];
    readonly singwebengine_get_total_reads: (a: number) => number;
    readonly singwebengine_get_mapped_reads: (a: number) => number;
    readonly singwebengine_get_mapped_bases: (a: number) => number;
    readonly singwebengine_get_genome_size: (a: number) => number;
    readonly singwebengine_get_avg_coverage: (a: number) => number;
    readonly singwebengine_map_reads_chunk: (a: number, b: number, c: number) => [number, number, number];
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __externref_table_dealloc: (a: number) => void;
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
