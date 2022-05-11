SDL2 skeleton for lab assignments 1–3 of the KTH course DH2323,
Computer Graphics and Interaction (and also SU DA3001, Datalogi II).

Build with:

    $ cmake .
    $ make

If successful, you should be able to run `./DgiLab` and see a
window completely filled with blue pixels.

# Requirements
You have to use a compiler supporting the C++11 standard. You also
need SDL2, GLM and CMake installed (see below).


### OS X
If you're using homebrew, I think running

    $ brew install sdl2 sdl2_gfx glm cmake

should be enough to get you going but I'm not sure as I'm also
using LLVM with Clang:

    $ brew install --with-clang llvm


### Linux
Markus Lidström reported that he had to switch out `GLM` for `OpenGL`
in `CMakeLists.txt` in order to make it run on his distribution of
Linux.
