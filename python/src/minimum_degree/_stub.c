/* The C library (../c/md.c) is built as an extension module so that pip
 * compiles it with the package; the wrapper loads it with ctypes, and never
 * imports it. */
#ifdef _WIN32
__declspec(dllexport)
#endif
void *PyInit__libmd(void) {
    return 0;
}
