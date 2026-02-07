import init, { SingWebEngine } from './pkg/sing_browser.js';

let engine = null;

self.onmessage = async (e) => {
    const { type, payload, id } = e.data;

    try {
        let result;
        switch (type) {
            case 'INIT':
                
                await init(new URL('./pkg/sing_browser_bg.wasm', import.meta.url));
                engine = new SingWebEngine();
                result = { success: true };
                break;
            case 'PROCESS_REF':
                if (!engine) throw new Error("Engine not initialized");
                engine.process_reference_chunk(payload.chunk, payload.chunkId);
                result = { success: true };
                break;
            case 'EXPORT_INDEX':
                if (!engine) throw new Error("Engine not initialized");
                const idxData = engine.export_index();
                result = { indexData: idxData, transfer: [idxData.buffer] };
                break;
            case 'IMPORT_INDEX':
                if (!engine) throw new Error("Engine not initialized");
                engine.import_index(payload.indexData);
                result = { success: true };
                break;
            case 'MAP_CHUNK':
                if (!engine) throw new Error("Engine not initialized");
                const { r1, r2, format, sort, writeHeader, filterMask } = payload;
                const mappingResult = engine.run_mapping_chunk(r1, r2, format, sort, writeHeader, filterMask);

                const bamBytes = mappingResult.bam_data;
                const mappedMask = mappingResult.mapped_mask;
                
                const transfer = [bamBytes.buffer, mappedMask.buffer];
                
                const genomeSize = engine.get_genome_size();
                
                result = { 
                    bam: bamBytes,
                    mappedMask,
                    stats: { 
                        total: mappingResult.total_reads, 
                        mapped: mappingResult.mapped_reads, 
                        mappedBases: mappingResult.mapped_bases, 
                        genomeSize 
                    },
                    transfer 
                };
                break;
            default:
                throw new Error(`Unknown message type: ${type}`);
        }
        
        
        if (result && result.transfer) {
             const transferList = result.transfer;
             delete result.transfer;
             self.postMessage({ type, id, result }, transferList);
        } else {
             self.postMessage({ type, id, result });
        }

    } catch (error) {
        self.postMessage({ type, id, error: error.message || error.toString() });
    }
};