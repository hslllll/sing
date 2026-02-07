// Export device detection for external use
export { detectDevice };

class WorkerClient {
    constructor() {
        this.worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
        this.promises = new Map();
        this.nextId = 0;
        this.busy = false;
        
        this.worker.onmessage = (e) => {
            const { type, id, result, error } = e.data;
            const p = this.promises.get(id);
            if (p) {
                if (error) p.reject(new Error(error));
                else p.resolve(result);
                this.promises.delete(id);
            }
        };
    }

    call(type, payload, transfer = []) {
        return new Promise((resolve, reject) => {
            const id = this.nextId++;
            this.promises.set(id, { resolve, reject });
            this.worker.postMessage({ type, payload, id }, transfer);
        });
    }
    
    init() { return this.call('INIT'); }
    
    processRef(chunk, chunkId) { 
        return this.call('PROCESS_REF', { chunk, chunkId }); 
    }

    exportIndex() {
        return this.call('EXPORT_INDEX', {});
    }

    importIndex(indexData) {
        return this.call('IMPORT_INDEX', { indexData });
    }
    
    mapChunk(r1, r2, format, sort, writeHeader, filterMask) { 
        this.busy = true;
        const transfer = [r1.buffer];
        if (r2) transfer.push(r2.buffer);
        if (filterMask) transfer.push(filterMask.buffer);
        return this.call('MAP_CHUNK', { r1, r2, format, sort, writeHeader, filterMask }, transfer)
            .finally(() => { this.busy = false; });
    }
    
    terminate() { this.worker.terminate(); }
}


// Amortized O(1) append buffer — eliminates O(n²) copy-on-append
class GrowableBuffer {
    constructor(initialCapacity = 256 * 1024) {
        this._buf = new Uint8Array(initialCapacity);
        this.length = 0;
    }
    _ensureCapacity(needed) {
        if (needed <= this._buf.length) return;
        let cap = this._buf.length;
        while (cap < needed) cap *= 2;
        const nb = new Uint8Array(cap);
        if (this.length > 0) nb.set(this._buf.subarray(0, this.length));
        this._buf = nb;
    }
    append(data) {
        this._ensureCapacity(this.length + data.length);
        this._buf.set(data, this.length);
        this.length += data.length;
    }
    consume(n) {
        const out = this._buf.slice(0, n);
        const rem = this.length - n;
        if (rem > 0) this._buf.copyWithin(0, n, this.length);
        this.length = rem;
        return out;
    }
    view() { return this._buf.subarray(0, this.length); }
}


function getFileStream(file) {
  const rawStream = file.stream();
  if (file.name.endsWith('.gz')) {
    if (typeof DecompressionStream !== 'undefined') {
        try {
            const decompressor = new DecompressionStream('gzip');
            return rawStream.pipeThrough(decompressor);
        } catch (e) {
            console.warn("DecompressionStream init failed, trying pako fallback", e);
        }
    }
    
    
    if (typeof pako !== 'undefined') {
        const inflator = new pako.Inflate();
        let buffer = new Uint8Array(0);
        let closed = false;

        return new ReadableStream({
            async start(controller) {
                const reader = rawStream.getReader();
                
                inflator.onData = (chunk) => {
                    controller.enqueue(chunk);
                };
                
                while (!closed) {
                    const { value, done } = await reader.read();
                    if (done) {
                        try {
                            inflator.push(new Uint8Array(0), true); 
                        } catch(err) { console.error("Inflate flush error", err); }
                        controller.close();
                        closed = true;
                        break;
                    }
                    if (value) {
                         try { 
                             inflator.push(value, false); 
                         } catch(err) {
                             controller.error(err);
                             closed = true;
                         }
                    }
                }
            }
        });
    }
    
    console.warn("No decompression method available. Returning raw stream (will likely fail for GZIP).");
    return rawStream;
  }
  return rawStream;
}



async function countReads(file, signal) {
    const reader = getFileStream(file).getReader();
    let newlineCount = 0;
    let bytesRead = 0;

    while (true) {
        if (signal?.aborted) throw new DOMException('aborted', 'AbortError');
        const { value, done } = await reader.read();
        if (value) {
            bytesRead += value.length;
            for (let i = 0; i < value.length; i++) {
                if (value[i] === 10) newlineCount++;
            }
            console.log(`Read ${bytesRead} bytes, ${Math.floor(newlineCount / 4)} reads counted so far...`);
        }
        if (done) break;
    }

    return Math.floor(newlineCount / 4);
}


function findSplitIndex(buffer, isFasta, capBytes = Infinity) {
    const len = buffer.length;

    if (len === 0) return -1;

    if (isFasta) {
        let lastHeaderStart = -1;
        
        for (let i = 0; i < len; i++) {
            if (buffer[i] === 10) { 
                if (i + 1 <= capBytes) {
                    
                    if (i + 1 < len && buffer[i+1] === 62) {
                         lastHeaderStart = i + 1;
                    }
                }
            }
            if (i + 1 > capBytes && lastHeaderStart !== -1) break;
        }

        
        if (lastHeaderStart !== -1) return lastHeaderStart;

        
        
        
        let lastNewline = -1;
        for (let i = 0; i < len; i++) {
            if (buffer[i] === 10) {
                if (i + 1 <= capBytes) lastNewline = i + 1;
            }
            if (i + 1 > capBytes && lastNewline !== -1) break;
        }
        return lastNewline;
    }

    
    let newlineCount = 0;
    let lastGood = -1;
    for (let i = 0; i < len; i++) {
        if (buffer[i] === 10) {
            newlineCount++;
            if (i + 1 <= capBytes && (newlineCount % 4 === 0)) {
                lastGood = i + 1;
            }
        }

        if (i + 1 > capBytes && lastGood !== -1) break;
    }

    return lastGood;
}

function recordsInFastq(buffer) {
        let newlineCount = 0;
        for (let i = 0; i < buffer.length; i++) {
                if (buffer[i] === 10) newlineCount++;
        }
        if (newlineCount === 0) return 0;
        if (newlineCount % 4 !== 0) return -1;
        return newlineCount / 4;
}

async function streamReadsPaired(file1, file2, chunkSize, onChunk) {
  const reader1 = getFileStream(file1).getReader();
  const reader2 = getFileStream(file2).getReader();
  
  const gbuf1 = new GrowableBuffer(chunkSize + 256 * 1024);
  const gbuf2 = new GrowableBuffer(chunkSize + 256 * 1024);
  
  let chunkId = 0;
  
  while (true) {
    
    const [r1, r2] = await Promise.all([reader1.read(), reader2.read()]);
    
    if (r1.value) gbuf1.append(r1.value);
    if (r2.value) gbuf2.append(r2.value);

    const done1 = r1.done;
    const done2 = r2.done;
    
    
    if (gbuf1.length >= chunkSize || gbuf2.length >= chunkSize || (done1 && done2)) {
        
        let split1 = findSplitIndex(gbuf1.view(), false, chunkSize);
        if (split1 === -1 && gbuf1.length > chunkSize * 4) {
            split1 = findSplitIndex(gbuf1.view(), false, Infinity);
        }
        
        if (split1 > 0) {
            const view1 = gbuf1.view();
            let records = 0;
            for(let i=0; i<split1; i++) if (view1[i]===10) records++;
            records = records / 4;
            
            let split2 = -1;
            let currentRecLines = 0;
            const view2 = gbuf2.view();
            for(let i=0; i<view2.length; i++) {
                if (view2[i]===10) {
                    currentRecLines++;
                    if (currentRecLines === records * 4) {
                        split2 = i + 1;
                        break;
                    }
                }
            }
            
            if (split2 === -1 && !done2) {
                if (done1 && r1.value === undefined) { 
                     if (r2.done) break; 
                     continue;
                }
                continue; 
            }
            
            const limit2 = (split2 !== -1) ? split2 : gbuf2.length;
            
            const chunk1 = gbuf1.consume(split1);
            const chunk2 = gbuf2.consume(limit2);
            
            await onChunk(chunk1, chunk2, chunkId, records);
            chunkId++;
            
            if (done1 && done2 && gbuf1.length === 0 && gbuf2.length === 0) break;
            
        } else {
             if (done1 && done2) {
                 const rec1 = recordsInFastq(gbuf1.view());
                 const rec2 = recordsInFastq(gbuf2.view());
                 if (rec1 > 0 && rec1 === rec2) {
                     await onChunk(gbuf1.consume(gbuf1.length), gbuf2.consume(gbuf2.length), chunkId, rec1);
                     chunkId++;
                 }
                 break; 
             }
        }
    } else {
        if (done1 && done2) break;
    }
  }
}

async function streamReadsSingle(file, chunkSize, onChunk) {
  const reader = getFileStream(file).getReader();
  const gbuf = new GrowableBuffer(chunkSize + 256 * 1024);
  let chunkId = 0;

  while(true) {
      const {value, done} = await reader.read();
      if (value) gbuf.append(value);
      
      if (gbuf.length >= chunkSize || done) {
          let split = findSplitIndex(gbuf.view(), false, chunkSize);
          if (split === -1 && gbuf.length > chunkSize * 4) {
              split = findSplitIndex(gbuf.view(), false, Infinity);
          }
          if (split > 0) {
              const view = gbuf.view();
              let records = 0;
              for (let i = 0; i < split; i++) if (view[i] === 10) records++;
              records = records / 4;

              const chunk = gbuf.consume(split);
              await onChunk(chunk, chunkId, records);
              chunkId++;
          }
          if (done && gbuf.length === 0) break;
          
          if (done && split === -1 && gbuf.length > 0) {
              const recs = recordsInFastq(gbuf.view());
              if (recs > 0) {
                  await onChunk(gbuf.consume(gbuf.length), chunkId, recs);
              }
              break;
          }
      }
      if (done && gbuf.length === 0) break;
  }
}


function concatPartsRange(parts, endPartIdx, endOffset) {
    let size = 0;
    for (let i = 0; i < endPartIdx; i++) size += parts[i].length;
    size += endOffset;
    
    const ret = new Uint8Array(size);
    let offset = 0;
    for (let i = 0; i < endPartIdx; i++) {
        ret.set(parts[i], offset);
        offset += parts[i].length;
    }
    if (endOffset > 0) {
        ret.set(parts[endPartIdx].subarray(0, endOffset), offset);
    }
    return ret;
}

function concatParts(parts, totalLen) {
    const ret = new Uint8Array(totalLen);
    let offset = 0;
    for(const p of parts) {
        ret.set(p, offset);
        offset += p.length;
    }
    return ret;
}


function preserveRemainder(parts, startPartIdx, startOffset) {
    const remains = [];
    
    if (startOffset < parts[startPartIdx].length) {
        remains.push(parts[startPartIdx].slice(startOffset)); 
        
    }
    
    for (let i = startPartIdx + 1; i < parts.length; i++) {
        remains.push(parts[i]);
    }
    return remains;
}


function findLastSplit(parts) {
    
    let absoluteIndex = 0; 
    
    
    
    if (parts.length === 0) return null;

    for (let i = parts.length - 1; i >= 0; i--) {
        const p = parts[i];
        for (let j = p.length - 1; j >= 0; j--) {
            if (p[j] === 62) { 
                
                let isNewline = false;
                if (j > 0) {
                    isNewline = (p[j-1] === 10);
                } else {
                    
                    if (i > 0) {
                        const prev = parts[i-1];
                        if (prev.length > 0) {
                            isNewline = (prev[prev.length - 1] === 10);
                        }
                    } else {
                        
                        
                        
                        
                        continue;
                    }
                }
                
                if (isNewline) {
                    return { partIndex: i, offsetInPart: j };
                }
            }
        }
    }
    return null;
}

async function streamReference(file, chunkSize, onChunk) {
    if (file.name.endsWith('.gz')) {
        throw new Error("GZIP args not supported. Uncompressed FASTA only.");
    }

    const reader = file.stream().getReader();
    let parts = []; 
    let partsTotalLength = 0;
    let chunkId = 0;
    
    while(true) {
        const {value, done} = await reader.read();
        if (value) {
            parts.push(value);
            partsTotalLength += value.length;
        }
        
        
        
        
        
        
        if (partsTotalLength >= chunkSize || done) {
             const split = findLastSplit(parts);
             
             if (split) {
                 
                 
                 const chunk = concatPartsRange(parts, split.partIndex, split.offsetInPart);
                 if (chunk.length > 0) {
                     await onChunk(chunk, chunkId++);
                 }
                 
                 
                 parts = preserveRemainder(parts, split.partIndex, split.offsetInPart);
                 partsTotalLength = parts.reduce((acc, p) => acc + p.length, 0);
                 
                 
                 
                 
                 
                 
                 
                 
                 
             } else {
                 if (done) {
                     
                     if (partsTotalLength > 0) {
                         const chunk = concatParts(parts, partsTotalLength);
                         await onChunk(chunk, chunkId++);
                     }
                     break; 
                 }
                 
                 
             }
        }
        
        if (done) break;
    }
}



const READ_CHUNK_SIZE = 100 * 1024 * 1024;


// Device detection and optimization
function detectDevice() {
    const mobile = /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent);
    const memory = navigator.deviceMemory || 8; // GB, fallback to 8
    const cores = navigator.hardwareConcurrency || 4;
    
    let profile = 'desktop';
    let recommendedChunkSize = 100 * 1024 * 1024; // 100MB
    let recommendedWorkers = cores;
    
    if (mobile) {
        if (memory <= 2) {
            profile = 'mobile-low';
            recommendedChunkSize = 10 * 1024 * 1024; // 10MB
            recommendedWorkers = Math.max(1, Math.floor(cores / 2));
        } else if (memory <= 4) {
            profile = 'mobile-mid';
            recommendedChunkSize = 25 * 1024 * 1024; // 25MB
            recommendedWorkers = Math.max(2, Math.floor(cores * 0.75));
        } else {
            profile = 'mobile-high';
            recommendedChunkSize = 50 * 1024 * 1024; // 50MB
            recommendedWorkers = cores;
        }
    }
    
    return { mobile, memory, cores, profile, recommendedChunkSize, recommendedWorkers };
}

// Performance metrics tracker
class PerformanceTracker {
    constructor() {
        this.startTime = Date.now();
        this.lastUpdate = Date.now();
        this.processedReads = 0;
        this.processedBytes = 0;
        this.readsPerSecond = 0;
        this.mbPerSecond = 0;
        // Self-tracked memory: accumulate known allocation sizes
        this.allocatedBytes = 0;
        this.peakAllocatedBytes = 0;
        // Chrome performance.memory (if supported)
        this._hasMemoryAPI = typeof performance !== 'undefined' && !!performance.memory;
        this.peakHeap = 0;
    }
    
    update(reads, bytes) {
        this.processedReads += reads;
        this.processedBytes += bytes;
        
        const now = Date.now();
        const elapsed = (now - this.lastUpdate) / 1000;
        
        if (elapsed > 0) {
            this.readsPerSecond = Math.round(reads / elapsed);
            this.mbPerSecond = (bytes / elapsed / (1024 * 1024)).toFixed(2);
        }
        
        this.lastUpdate = now;
        
        if (this._hasMemoryAPI) {
            this.peakHeap = Math.max(this.peakHeap, performance.memory.usedJSHeapSize);
        }
    }
    
    // Call when large allocations are created/released
    trackAlloc(bytes) {
        this.allocatedBytes += bytes;
        this.peakAllocatedBytes = Math.max(this.peakAllocatedBytes, this.allocatedBytes);
    }
    trackFree(bytes) {
        this.allocatedBytes = Math.max(0, this.allocatedBytes - bytes);
    }
    
    getStats() {
        const elapsed = (Date.now() - this.startTime) / 1000;
        const avgReadsPerSec = elapsed > 0 ? Math.round(this.processedReads / elapsed) : 0;
        const avgMbPerSec = elapsed > 0 ? (this.processedBytes / elapsed / (1024 * 1024)).toFixed(2) : 0;
        
        let currentMemoryMB, peakMemoryMB;
        if (this._hasMemoryAPI) {
            currentMemoryMB = (performance.memory.usedJSHeapSize / (1024 * 1024)).toFixed(1);
            peakMemoryMB = (this.peakHeap / (1024 * 1024)).toFixed(1);
        } else {
            // Fallback: report tracked buffer sizes
            currentMemoryMB = (this.allocatedBytes / (1024 * 1024)).toFixed(1);
            peakMemoryMB = (this.peakAllocatedBytes / (1024 * 1024)).toFixed(1);
        }
        
        return {
            elapsed: elapsed.toFixed(1),
            processedReads: this.processedReads,
            processedMB: (this.processedBytes / (1024 * 1024)).toFixed(2),
            avgReadsPerSec,
            avgMbPerSec,
            currentReadsPerSec: this.readsPerSecond,
            currentMbPerSec: this.mbPerSecond,
            peakMemoryMB,
            currentMemoryMB,
            memorySource: this._hasMemoryAPI ? 'heap' : 'tracked'
        };
    }
}


export async function partitionedWorkflow(referenceFile, read1File, read2File, sortOutput = false, outputFormat = 'bam', chunkSize = READ_CHUNK_SIZE, logger = console.log, onProgress = null, onMetrics = null) {
    const abortCtrl = new AbortController();
    let totalReads = 0;
    const formatTime = ms => !isFinite(ms) || ms < 0 ? "--:--" : `${Math.floor(ms / 60000)}:${(Math.floor(ms / 1000) % 60).toString().padStart(2, '0')}`;

    // Detect device and optimize settings
    const device = detectDevice();
    logger(`Device: ${device.profile} | Memory: ${device.memory}GB | Cores: ${device.cores}`);
    
    // Auto-adjust chunk size if default is being used
    if (chunkSize === READ_CHUNK_SIZE && device.mobile) {
        chunkSize = device.recommendedChunkSize;
        logger(`Auto-adjusted chunk size to ${(chunkSize / (1024 * 1024)).toFixed(0)}MB for mobile`);
    }
    
    // Initialize performance tracker
    const perfTracker = new PerformanceTracker();

    // Quick file-size estimate — avoids reading entire files before alignment
    const estimateReads = (f) => {
        const eff = f.name.endsWith('.gz') ? f.size * 3.5 : f.size;
        return Math.floor(eff / 400); // ~400 bytes per FASTQ record
    };
    totalReads = estimateReads(read1File);
    if (read2File) totalReads = Math.min(totalReads, estimateReads(read2File));
    logger(`Estimated reads: ~${totalReads.toLocaleString()}`);

    const concurrency = device.mobile ? device.recommendedWorkers : (navigator.hardwareConcurrency || 4);
    logger(`Initializing protocol with ${concurrency} workers...`);
    const workers = Array.from({ length: concurrency }, () => new WorkerClient());
    await Promise.all(workers.map(w => w.init()));

    const runStats = { mapped: 0, mappedBases: 0, genomeSize: 0, total: 0 };
    const finalTotalReads = totalReads * (read2File ? 2 : 1);
    
    const mappedMaskSize = Math.ceil((finalTotalReads || 1) / 8);
    const globalMappedMask = new Uint8Array(mappedMaskSize);
    perfTracker.trackAlloc(mappedMaskSize);
    
    const bamParts = [];
    let processedReads = 0, processedBytes = 0;
    const startTime = Date.now();
    const refChunkSize = 50 * 1024 * 1024;
    const estimatedRefChunks = Math.max(1, Math.ceil(referenceFile.size / refChunkSize));
    const totalProgressWork = Math.max(finalTotalReads, 1) * estimatedRefChunks;

    const updateProgress = (newReads, newBytes) => {
        processedReads += newReads;
        processedBytes += newBytes;
        
        perfTracker.update(newReads, newBytes);
        const stats = perfTracker.getStats();
        
        if (onMetrics) {
            onMetrics(stats);
        }
        
        const elapsed = Date.now() - startTime;
        return `Elapsed: ${formatTime(elapsed)} | ${stats.currentReadsPerSec.toLocaleString()} reads/s | ${stats.currentMbPerSec} MB/s | Mem: ${stats.currentMemoryMB} MB`; 
    };

    try {
        await streamReference(referenceFile, refChunkSize, async (refChunk, refId) => {
            logger(`RefChunk ${refId}: ${(refChunk.length/1e6).toFixed(2)} MB`);
            perfTracker.trackAlloc(refChunk.byteLength);
            
            await workers[0].processRef(refChunk, refId);
            if (workers.length > 1) {
                const { indexData } = await workers[0].exportIndex();
                perfTracker.trackAlloc(indexData.byteLength);
                await Promise.all(workers.slice(1).map(w => w.importIndex(indexData)));
                perfTracker.trackFree(indexData.byteLength);
            }
            perfTracker.trackFree(refChunk.byteLength);
            
            const activePromises = new Set();
            const chunkParts = []; 
            let passGenomeSize = 0;
            let readChunkCount = 0;
            let currentReadOffset = 0;

            const handleResult = (result, chunkId, chunkStartReadIdx) => {
                const { bam, stats, mappedMask } = result;
                runStats.mapped += stats.mapped; 
                runStats.mappedBases += stats.mappedBases;
                if (refId === 0) runStats.total += stats.total;
                passGenomeSize = Math.max(passGenomeSize, stats.genomeSize);
                
                const startBit = chunkStartReadIdx;
                for(let i=0; i<mappedMask.length; i++) {
                    const byteVal = mappedMask[i];
                    if (byteVal === 0) continue;
                    
                    for(let b=0; b<8; b++) {
                        if ((byteVal & (1 << b)) !== 0) {
                            const globalBit = startBit + i*8 + b;
                            const globalByte = Math.floor(globalBit / 8);
                            const globalBitInByte = globalBit % 8;
                            if (globalByte < globalMappedMask.length) {
                                globalMappedMask[globalByte] |= (1 << globalBitInByte);
                            }
                        }
                    }
                }
                

                if (bam.length) {
                    chunkParts[chunkId] = new Uint8Array(bam);
                }
            };
            
            const processChunk = async (r1, r2, cId, recs) => {
                readChunkCount++;
                const nReads = recs * (r2 ? 2 : 1);
                const thisChunkStartIdx = currentReadOffset;
                currentReadOffset += nReads;

                while (activePromises.size >= workers.length) await Promise.race(activePromises);
                
                const w = workers.find(x => !x.busy);
                const writeHead = (readChunkCount === 1 && refId === 0);
                
                const chunkMaskSize = Math.ceil(nReads / 8);
                const chunkFilterMask = new Uint8Array(chunkMaskSize);
                
                for(let i=0; i<nReads; i++) {
                    const globalBit = thisChunkStartIdx + i;
                    const globalByte = Math.floor(globalBit / 8);
                    const globalBitInByte = globalBit % 8;
                    const isSet = (globalMappedMask[globalByte] & (1 << globalBitInByte)) !== 0;
                    
                    if (isSet) {
                        const localByte = Math.floor(i / 8);
                        const localBit = i % 8;
                        chunkFilterMask[localByte] |= (1 << localBit);
                    }
                }
                
                // Capture byte sizes BEFORE transfer (buffers get detached)
                const r1Bytes = r1.byteLength;
                const r2Bytes = r2 ? r2.byteLength : 0;
                const chunkBytes = r1Bytes + r2Bytes;
                perfTracker.trackAlloc(chunkBytes);
                
                const p = w.mapChunk(r1, r2, outputFormat, sortOutput, writeHead, chunkFilterMask)
                    .then(res => {
                        perfTracker.trackFree(chunkBytes);
                        handleResult(res, cId, thisChunkStartIdx);
                        perfTracker.trackAlloc(res.bam.byteLength);
                        const msg = updateProgress(nReads, r1Bytes + r2Bytes);
                        logger(msg);
                        if(onProgress) onProgress(((processedReads / totalProgressWork)*100).toFixed(1));
                    });
                
                activePromises.add(p);
                p.finally(() => activePromises.delete(p));
            };

            if (read2File) await streamReadsPaired(read1File, read2File, chunkSize, processChunk);
            else await streamReadsSingle(read1File, chunkSize, (r1, cId, recs) => processChunk(r1, undefined, cId, recs));
            
            await Promise.all(activePromises);
            runStats.genomeSize += passGenomeSize; 
            
            const globalCov = runStats.genomeSize ? (runStats.mappedBases / runStats.genomeSize).toFixed(2) : 0;
            logger(`RefChunk ${refId} done. Global Mapped: ${runStats.mapped}, Genome: ${runStats.genomeSize}, Avg Coverage: ${globalCov}x`);

            chunkParts.forEach(p => { if (p) bamParts.push(p); });
        });

        logger("Merging...");
        if (outputFormat === 'bam') {
            bamParts.push(new Uint8Array([0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00])); 
        }

        const totalTime = Date.now() - startTime;
        const finalStats = perfTracker.getStats();
        logger(`Done in ${formatTime(totalTime)}.`);
        logger(`Performance: ${finalStats.avgReadsPerSec.toLocaleString()} avg reads/s | ${finalStats.avgMbPerSec} avg MB/s`);
        logger(`Memory: Peak ${finalStats.peakMemoryMB} MB | Current ${finalStats.currentMemoryMB} MB`);
        logger(`Processed: ${finalStats.processedReads.toLocaleString()} reads | ${finalStats.processedMB} MB`);
        
        if (onProgress) onProgress(100);
        if (onMetrics) onMetrics(finalStats);

        return { 
            blob: new Blob(bamParts, { type: 'application/octet-stream' }), 
            stats: { 
                total: runStats.total || finalTotalReads, 
                mapped: runStats.mapped, 
                mappedBases: runStats.mappedBases, 
                genomeSize: runStats.genomeSize, 
                avgCoverage: runStats.genomeSize ? runStats.mappedBases / runStats.genomeSize : 0 
            } 
        };
    } finally {
        workers.forEach(w => w.terminate());
    }
}