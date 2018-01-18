# Default bindings
bindgen --whitelist-function "ovr.*" --whitelist-function "vrapi.*" --whitelist-type "ovr.*" --whitelist-var "VRAPI.*" -o src/bindings.rs VrApi/wrapper.h -- -std=c99 -I/usr/include/clang/3.9/include

# Android bindings
ANDROID_INCLUDES="$ANDROID_NDK/platforms/android-18/arch-arm/usr/include"
bindgen --whitelist-function "ovr.*" --whitelist-function "vrapi.*" --whitelist-type "ovr.*" --whitelist-var "VRAPI.*" -o src/bindings_android.rs VrApi/wrapper.h -- -std=c99 -D__ANDROID__ -DANDROID -I$ANDROID_INCLUDES -I/usr/include/clang/3.9/include
