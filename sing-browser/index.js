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
    logger("Counting reads for progress tracking (skip with ESC)...");
    const abortCtrl = new AbortController();
    let aborted = false;
    const onKey = (e) => {
        if (e.key === 'Escape') {
            aborted = true;
            abortCtrl.abort();
            logger('Read counting skipped by user; using byte-based progress.');
            window.removeEventListener('keydown', onKey);
        }
    };
    window.addEventListener('keydown', onKey);

    let totalReadsR1 = 0;
    let totalReadsR2 = 0;
    try {
        totalReadsR1 = await countReads(read1File, abortCtrl.signal);
        logger(`Total reads in R1: ${totalReadsR1}`);
        if (read2File) {
            totalReadsR2 = await countReads(read2File, abortCtrl.signal);
            logger(`Total reads in R2: ${totalReadsR2}`);
            if (totalReadsR1 !== totalReadsR2) {
                logger(`Warning: read count mismatch between pairs (R1=${totalReadsR1}, R2=${totalReadsR2}). Progress will use the smaller value.`);
            }
        }
    } catch (e) {
        if (aborted) {
            totalReadsR1 = 0;
            totalReadsR2 = 0;
        } else {
            logger(`Read counting failed: ${e.message || e}`);
            totalReadsR1 = 0;
            totalReadsR2 = 0;
        }
    } finally {
        window.removeEventListener('keydown', onKey);
    }

    const baseReads = read2File ? Math.min(totalReadsR1, totalReadsR2) : totalReadsR1;
    const readsPerIteration = baseReads * (read2File ? 2 : 1);
    if (baseReads > 0) {
        logger(`Detected approximately ${readsPerIteration} reads per reference chunk iteration.`);
    }

    
    const concurrency = navigator.hardwareConcurrency || 4;
    logger(`Initializing Worker Pool with ${concurrency} workers...`);
    
    const workers = [];
    for(let i=0; i<concurrency; i++) {
        workers.push(new WorkerClient());
    }

    const workerStats = new Map();

    try {
        await Promise.all(workers.map(w => w.init()));
        logger("Workers Initialized.");
    
        const bamParts = [];
              let latestStats = null;
        logger("Calculating total work for progress tracking...");
        
          const refChunkSizeBytes = 100 * 1024 * 1024; 
        
        const estimatedRefChunks = Math.max(1, Math.ceil(referenceFile.size / refChunkSizeBytes));
          const totalWorkReads = readsPerIteration * estimatedRefChunks;
          const totalReadBytes = read1File.size + (read2File ? read2File.size : 0);
          const totalWorkBytes = totalReadBytes * estimatedRefChunks;
          let processedReads = 0;
          let processedBytes = 0;
      
          if (totalWorkReads === 0 && onProgress) {
                  logger("Read counting failed or yielded zero; falling back to byte-based progress.");
          }
        
          const startTime = Date.now();
      
          if (onProgress) {
                  onProgress(0);
          }
      
        function formatTime(ms) {
            if (!isFinite(ms) || ms < 0) return "--:--";
            const seconds = Math.floor(ms / 1000);
            const m = Math.floor(seconds / 60);
            const s = seconds % 60;
            return `${m}:${s.toString().padStart(2, '0')}`;
        }
      
        async function updateProgress(newlyProcessedReads, newlyProcessedBytes) {
            processedReads += newlyProcessedReads;
            processedBytes += newlyProcessedBytes;
            const elapsed = Date.now() - startTime;
      
            let fraction = 0;
            if (totalWorkReads > 0 && processedReads > 0) {
                fraction = processedReads / totalWorkReads;
            } else if (totalWorkBytes > 0 && processedBytes > 0) {
                fraction = processedBytes / totalWorkBytes;
            }
      
            if (fraction > 0) {
                const estimatedTotalTime = elapsed / fraction;
                const remaining = estimatedTotalTime - elapsed;
                const percentage = Math.min(100, fraction * 100);
      
                if (onProgress) {
                    onProgress(parseFloat(percentage.toFixed(1)));
                    
                    await flushUi();
                }
      
                return `[${percentage.toFixed(1)}%] Elapsed: ${formatTime(elapsed)}, ETA: ${formatTime(remaining)}`;
            }
            return ``;
        }
      
        
        
        await streamReference(referenceFile, refChunkSizeBytes, async (refChunk, refId) => {
            logger(`Processing Reference Chunk ${refId}, size: ${(refChunk.length / 1024 / 1024).toFixed(2)} MB`);
            
            
            if (workers.length > 1) {
                logger(`Building index in sub-worker 0...`);
                await workers[0].processRef(refChunk, refId);
                
                logger(`Broadcasting index to ${workers.length - 1} workers...`);
                const { indexData } = await workers[0].exportIndex();
                
                
                
                
                
                
                
                const promises = [];
                for (let i = 1; i < workers.length; i++) {
                    promises.push(workers[i].importIndex(indexData));
                }
                await Promise.all(promises);
            } else {
                await workers[0].processRef(refChunk, refId);
            }
            
            const activePromises = new Set();
            let readChunkCount = 0;

            const chunkParts = [];
            let maxChunkId = -1;

            const handleResult = (result, chunkId, worker) => {
                const { bam, stats } = result;
                
                if (worker) {
                    workerStats.set(worker, stats);
                }

                let aggTotal = 0;
                let aggMapped = 0;
                let aggMappedBases = 0;
                let aggGenomeSize = stats.genomeSize; 
                
                for (const s of workerStats.values()) {
                    aggTotal += s.total;
                    aggMapped += s.mapped;
                    aggMappedBases += s.mappedBases;
                }

                let aggAvgCoverage = 0;
                if (aggGenomeSize > 0) {
                    aggAvgCoverage = aggMappedBases / aggGenomeSize;
                }

                latestStats = {
                    total: aggTotal,
                    mapped: aggMapped,
                    mappedBases: aggMappedBases,
                    genomeSize: aggGenomeSize,
                    avgCoverage: aggAvgCoverage
                };

                logger(`Processed Chunk ${chunkId}: ${latestStats.mapped}/${latestStats.total} reads mapped.`);

                
                
                

                if (bam.length > 0) {
                    let part = new Uint8Array(bam);
                    
                    
                    if (outputFormat === 'sam' && bamParts.length > 0) {
                        
                        
                        
                        let splitIndex = 0;
                        const len = part.length;
                        let inHeader = true;
                        
                        
                        if (len > 0 && part[0] === 64) {
                            
                            for(let i=0; i<len; i++) {
                               if (part[i] === 10) { 
                                   if (i + 1 < len) {
                                       if (part[i+1] !== 64) {
                                           
                                           splitIndex = i + 1;
                                           inHeader = false;
                                           break;
                                       }
                                   } else {
                                       
                                       splitIndex = len;
                                       inHeader = false;
                                   }
                               }
                            }
                            if (inHeader) splitIndex = len;
                        }
                        
                        if (splitIndex > 0) {
                            part = part.subarray(splitIndex);
                        }
                    }
                    
                    if (part.length > 0) {
                        chunkParts[chunkId] = part;
                        if (chunkId > maxChunkId) maxChunkId = chunkId;
                    }
                }
            };
            
            if (read2File) {
                logger(`Mapping Paired-End reads against Reference Chunk ${refId}...`);
                
                await streamReadsPaired(read1File, read2File, chunkSize, async (r1, r2, chunkId, records) => {
                     readChunkCount++;
                     const currentChunkReads = records * 2;
                     const currentChunkBytes = r1.byteLength + r2.byteLength;
                     
                     
                     while (activePromises.size >= workers.length) {
                         await Promise.race(activePromises);
                     }
                     
                     const worker = workers.find(w => !w.busy);
                     if (!worker) throw new Error("Worker pool logic error");
                     
                     const writeHeader = (readChunkCount === 1) && (refId === 0);
                     const promise = worker.mapChunk(r1, r2, outputFormat, sortOutput, writeHeader);
                     activePromises.add(promise);
                     
                     const p = promise.then(result => {
                         activePromises.delete(promise);
                         handleResult(result, chunkId, worker);
                         return updateProgress(currentChunkReads, currentChunkBytes);
                     }).then(timeLog => {
                         
                     });
                     
                });
            } else {
                logger(`Mapping Single-End reads against Reference Chunk ${refId}...`);
                
                await streamReadsSingle(read1File, chunkSize, async (r1, chunkId, records) => {
                     readChunkCount++;
                     const currentChunkReads = records;
                     const currentChunkBytes = r1.byteLength;
                     
                     while (activePromises.size >= workers.length) {
                         await Promise.race(activePromises);
                     }
                     
                     const worker = workers.find(w => !w.busy);
                     const writeHeader = (readChunkCount === 1) && (refId === 0);
                     const promise = worker.mapChunk(r1, undefined, outputFormat, sortOutput, writeHeader);
                     activePromises.add(promise);
                     
                     const p = promise.then(result => {
                         activePromises.delete(promise);
                         handleResult(result, chunkId, worker);
                         return updateProgress(currentChunkReads, currentChunkBytes);
                     }).then(timeLog => {
                         
                     });
                });
            }
            
            
            await Promise.all(activePromises);
            logger(`Completed mapping all reads against Reference Chunk ${refId}. Total Read Chunks: ${readChunkCount}`);

            
            for (let i = 0; i <= maxChunkId; i++) {
                const part = chunkParts[i];
                if (part && part.length > 0) {
                    bamParts.push(part);
                    chunkParts[i] = null; // release chunk memory early
                }
            }
        });
      
        logger("Merging BAM parts...");
        
        if (outputFormat === 'bam') {
            
            const BAM_EOF = new Uint8Array([
                0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 
                0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
            ]);
            bamParts.push(BAM_EOF);
        }

        const blob = new Blob(bamParts, { type: 'application/octet-stream' });
        
        const totalTime = Date.now() - startTime;
        logger(`Total Mapping Time: ${formatTime(totalTime)} ( ${(totalTime/1000).toFixed(2)} seconds)`);

        if (onProgress) {
            onProgress(100);
        }

        return { blob, stats: latestStats };
    } finally {
        
        workers.forEach(w => w.terminate());
    }
}