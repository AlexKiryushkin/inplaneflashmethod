
#pragma once

#ifdef IPFM_EXPORTS
#    define IPFM_API __declspec(dllexport)
#else
#    define IPFM_API __declspec(dllimport)
#endif
