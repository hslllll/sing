#[cfg(target_arch = "wasm32")]
mod wasm;

#[cfg(target_arch = "wasm32")]
pub use wasm::SingWebEngine;

#[cfg(not(target_arch = "wasm32"))]
pub struct SingWebEngine;

#[cfg(not(target_arch = "wasm32"))]
impl SingWebEngine {
    pub fn new() -> Self {
        SingWebEngine
    }

    pub fn process_reference_chunk(&mut self, _chunk_data: &[u8], _chunk_id: usize) -> Result<(), String> {
        Err("sing-web is only available when compiled for wasm32-unknown-unknown".to_string())
    }

    pub fn map_reads_chunk(&mut self, _reads_chunk: &[u8]) -> Result<Vec<u8>, String> {
        Err("sing-web is only available when compiled for wasm32-unknown-unknown".to_string())
    }
}
