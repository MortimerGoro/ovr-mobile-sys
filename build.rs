extern crate bindgen;

use std::env;
use std::path::{Path, PathBuf};
use bindgen::Builder;

fn main() {
    let target = env::var("TARGET").unwrap();
    let android = target.contains("android");
 
    // Generate gvr Rust bindings.
    let mut bindings = Builder::default()
                .no_unstable_rust()
                .header("VrApi/Include/VrApi.h")
                .header("VrApi/Include/VrApi_Input.h")
                .header("VrApi/Include/VrApi_Version.h")
                .header("VrApi/Include/VrApi_Helpers.h")
                .clang_arg("-std=c99");

    if android {
        // Add Android macros and sysroot to generate correct JNI methods in the API
        let ndk_path = env::var("ANDROID_NDK").ok().expect("Please set the ANDROID_NDK environment variable");
        let ndk_path = Path::new(&ndk_path);
        let (arch, platform) = if target.contains("x86") {
            ("arch-x86", "android-18")
        } else if target.contains("aarch64") {
            ("arch-arm64", "android-21")
        } else {
            ("arch-arm", "android-18")
        };
        let sys_root = ndk_path.join("platforms").join(platform).join(arch).join("usr").join("include");
        bindings = bindings.clang_arg("-D__ANDROID__")
                           .clang_arg(format!("-I{}", sys_root.to_str().unwrap()));
    }
            
    let bindings = bindings.generate()
                           .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings.write_to_file(out_path.join("bindings.rs"))
            .expect("Couldn't write bindings!");

    // Export shared libraries search path.
    if android {
        let abi = if target.contains("aarch64") {
            "arm64-v8a"
        } else {
            "armeabi-v7a"
        };
        
        println!("cargo:rustc-link-search={}/VrApi/Libs/Android/{}", env!("CARGO_MANIFEST_DIR"), abi);
        println!("cargo:rustc-link-lib=vrapi");
    }
}
