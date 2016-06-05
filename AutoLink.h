#ifndef AUTOLINK_H_
#define AUTOLINK_H_

#define NAME "MG"

#if defined(_DEBUG)
#	define TYPE "-d"
#else
#   define TYPE "-r"
#endif

#define LIB_NAME NAME TYPE ".lib"

#pragma comment(lib, LIB_NAME)

#endif