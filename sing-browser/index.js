import init, { SingWebEngine } from './pkg/sing_browser.js';


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
    
    mapChunk(r1, r2, format, sort, writeHeader) { 
        this.busy = true;
        const transfer = [r1.buffer];
        if (r2) transfer.push(r2.buffer);
        return this.call('MAP_CHUNK', { r1, r2, format, sort, writeHeader }, transfer)
            .finally(() => { this.busy = false; });
    }
    
    terminate() { this.worker.terminate(); }
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
  
  let buf1 = new Uint8Array(0);
  let buf2 = new Uint8Array(0);
  
  let chunkId = 0;
  
  while (true) {
    
    const [r1, r2] = await Promise.all([reader1.read(), reader2.read()]);
    
    
    if (r1.value) {
        const next = new Uint8Array(buf1.length + r1.value.length);
        next.set(buf1);
        next.set(r1.value, buf1.length);
        buf1 = next;
    }
    if (r2.value) {
        const next = new Uint8Array(buf2.length + r2.value.length);
        next.set(buf2);
        next.set(r2.value, buf2.length);
        buf2 = next;
    }

    const done1 = r1.done;
    const done2 = r2.done;
    
    
    if (buf1.length >= chunkSize || buf2.length >= chunkSize || (done1 && done2)) {
        
        
        
        
        let split1 = findSplitIndex(buf1, false, chunkSize);
        if (split1 === -1 && buf1.length > chunkSize * 4) {
            
            split1 = findSplitIndex(buf1, false, Infinity);
        }
        
        if (split1 > 0) {
            
            let records = 0;
            for(let i=0; i<split1; i++) if (buf1[i]===10) records++;
            records = records / 4;
            
            
            let split2 = -1;
            let currentRecLines = 0;
            for(let i=0; i<buf2.length; i++) {
                if (buf2[i]===10) {
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
            
            
            const limit2 = (split2 !== -1) ? split2 : buf2.length;
            
            const chunk1 = buf1.slice(0, split1);
            const chunk2 = buf2.slice(0, limit2);
            
            await onChunk(chunk1, chunk2, chunkId, records);
            chunkId++;
            
            buf1 = buf1.slice(split1);
            buf2 = buf2.slice(limit2);
            
            
            if (done1 && done2 && buf1.length === 0 && buf2.length === 0) break;
            
        } else {
             
             if (done1 && done2) {
                 const rec1 = recordsInFastq(buf1);
                 const rec2 = recordsInFastq(buf2);
                 if (rec1 > 0 && rec1 === rec2) {
                     await onChunk(buf1, buf2, chunkId, rec1);
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
  let buffer = new Uint8Array(0);
  let chunkId = 0;

  while(true) {
      const {value, done} = await reader.read();
      if (value) {
          const next = new Uint8Array(buffer.length + value.length);
          next.set(buffer);
          next.set(value, buffer.length);
          buffer = next;
      }
      
      if (buffer.length >= chunkSize || done) {
          let split = findSplitIndex(buffer, false, chunkSize);
          if (split === -1 && buffer.length > chunkSize * 4) {
              split = findSplitIndex(buffer, false, Infinity);
          }
          if (split > 0) {
              let records = 0;
              for (let i = 0; i < split; i++) if (buffer[i] === 10) records++;
              records = records / 4;

              const chunk = buffer.slice(0, split);
              await onChunk(chunk, chunkId, records);
              chunkId++;
              buffer = buffer.slice(split);
          }
          if (done && buffer.length === 0) break;
          
          if (done && split === -1 && buffer.length > 0) {
              const recs = recordsInFastq(buffer);
              if (recs > 0) {
                  await onChunk(buffer, chunkId, recs);
              }
              break;
          }
      }
      if (done && buffer.length === 0) break;
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



const CHUNK_SIZE = 1 * 1024 * 1024;


export async function partitionedWorkflow(referenceFile, read1File, read2File, sortOutput = false, outputFormat = 'bam', chunkSize = CHUNK_SIZE, logger = console.log, onProgress = null) {
    const abortCtrl = new AbortController();
    let totalReads = 0;
    const formatTime = ms => !isFinite(ms) || ms < 0 ? "--:--" : `${Math.floor(ms / 60000)}:${(Math.floor(ms / 1000) % 60).toString().padStart(2, '0')}`;

    try {
        logger("Counting reads...");
        const keyHandler = e => { if (e.key === 'Escape') { abortCtrl.abort(); logger('Count skipped.'); } };
        window.addEventListener('keydown', keyHandler);
        
        const r1Count = await countReads(read1File, abortCtrl.signal).catch(() => 0);
        let r2Count = r1Count;
        if (read2File) r2Count = await countReads(read2File, abortCtrl.signal).catch(() => 0);
        
        window.removeEventListener('keydown', keyHandler);
        totalReads = Math.min(r1Count, r2Count || r1Count);
        logger(`Approx reads: ${totalReads}`);
    } catch(e) { logger("Count failed."); }

    const concurrency = navigator.hardwareConcurrency || 4;
    logger(`Initializing protocol with ${concurrency} workers...`);
    const workers = Array.from({ length: concurrency }, () => new WorkerClient());
    await Promise.all(workers.map(w => w.init()));

    const runStats = { total: 0, mapped: 0, mappedBases: 0, genomeSize: 0 };
    const bamParts = [];
    let processedReads = 0, processedBytes = 0;
    const startTime = Date.now();
    const totalBytes = read1File.size + (read2File?.size || 0);

    const updateProgress = (newReads, newBytes) => {
        processedReads += newReads;
        processedBytes += newBytes;
        const elapsed = Date.now() - startTime;
        let p = 0;
        return `Elapsed: ${formatTime(elapsed)}`; 
    };

    try {
        const refChunkSize = 100 * 1024 * 1024;
        await streamReference(referenceFile, refChunkSize, async (refChunk, refId) => {
            logger(`RefChunk ${refId}: ${(refChunk.length/1e6).toFixed(2)} MB`);
            
            await workers[0].processRef(refChunk, refId);
            if (workers.length > 1) {
                const { indexData } = await workers[0].exportIndex();
                await Promise.all(workers.slice(1).map(w => w.importIndex(indexData)));
            }
            
            const activePromises = new Set();
            const chunkParts = []; 
            let passGenomeSize = 0;
            let readChunkCount = 0;

            const handleResult = (result, chunkId) => {
                const { bam, stats } = result;
                runStats.total += stats.total; 
                runStats.mapped += stats.mapped;
                runStats.mappedBases += stats.mappedBases;
                passGenomeSize = Math.max(passGenomeSize, stats.genomeSize);
                
                const cov = passGenomeSize ? (runStats.mappedBases / passGenomeSize).toFixed(2) : 0;
                logger(`Chunk ${chunkId}: ${stats.mapped}/${stats.total} mapped. Global Mapped: ${runStats.mapped}. Cov: ${cov}`);

                if (bam.length) {
                    chunkParts[chunkId] = new Uint8Array(bam);
                }
            };
            
            const processChunk = async (r1, r2, cId, recs) => {
                readChunkCount++;
                while (activePromises.size >= workers.length) await Promise.race(activePromises);
                
                const w = workers.find(x => !x.busy);
                const writeHead = (readChunkCount === 1 && refId === 0);
                
                const p = w.mapChunk(r1, r2, outputFormat, sortOutput, writeHead)
                    .then(res => {
                        handleResult(res, cId);
                        const msg = updateProgress(recs * (r2?2:1), r1.byteLength + (r2?r2.byteLength:0));
                        if(onProgress) onProgress(((processedReads / (totalReads || 1))*10).toFixed(1));
                    });
                
                activePromises.add(p);
                p.finally(() => activePromises.delete(p));
            };

            if (read2File) await streamReadsPaired(read1File, read2File, chunkSize, processChunk);
            else await streamReadsSingle(read1File, chunkSize, (r1, cId, recs) => processChunk(r1, undefined, cId, recs));
            
            await Promise.all(activePromises);
            runStats.genomeSize += passGenomeSize; 

            chunkParts.forEach(p => { if (p) bamParts.push(p); });
        });

        logger("Merging...");
        if (outputFormat === 'bam') {
            bamParts.push(new Uint8Array([0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00])); 
        }

        const totalTime = Date.now() - startTime;
        logger(`Done in ${formatTime(totalTime)}.`);
        if (onProgress) onProgress(100);

        return { 
            blob: new Blob(bamParts, { type: 'application/octet-stream' }), 
            stats: { ...runStats, avgCoverage: runStats.genomeSize ? runStats.mappedBases / runStats.genomeSize : 0 } 
        };
    } finally {
        workers.forEach(w => w.terminate());
    }
}