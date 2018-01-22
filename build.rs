use std::env;

fn main() {
    let target = env::var("TARGET").unwrap();
    let android = target.contains("android");

    // Export shared libraries search path.
    if android {
        let abi = if target.contains("aarch64") {
            "arm64-v8a"
        } else {
            "armeabi-v7a"
        };
        
        println!("cargo:rustc-link-search={}/VrApi/Libs/Android/{}/Release", env!("CARGO_MANIFEST_DIR"), abi);
        println!("cargo:rustc-link-lib=vrapi");
    }
}
