import init, { SingWebEngine } from './pkg/sing_web.js';

let engine = null;

self.onmessage = async (e) => {
    const { type, payload, id } = e.data;

    try {
        let result;
        switch (type) {
            case 'INIT':
                
                await init(new URL(`./pkg/sing_web_bg.wasm?t=${Date.now()}`, import.meta.url));
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
                const { r1, r2, format, sort, writeHeader } = payload;
                const mappingResult = engine.run_mapping_chunk(0, r1, r2, format, sort, writeHeader);

                const bamBytes = mappingResult.bam_data;
                const transfer = [bamBytes.buffer];

                
                const total = engine.get_total_reads();
                const mapped = engine.get_mapped_reads();
                const mappedBases = engine.get_mapped_bases();
                const genomeSize = engine.get_genome_size();
                const avgCoverage = engine.get_avg_coverage();
                
                result = { 
                    bam: bamBytes,
                    stats: { total, mapped, mappedBases, genomeSize, avgCoverage },
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