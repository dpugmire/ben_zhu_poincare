#ifndef CODEX_CXX2_VISKORESCOMPAT_H
#define CODEX_CXX2_VISKORESCOMPAT_H

#include <cstddef>

#if defined(CODEX_USE_VISKORES)
#include <viskores/Types.h>
#define CODEX_EXEC VISKORES_EXEC
#define CODEX_CONT VISKORES_CONT
using CodeXId = viskores::Id;
#else
#define CODEX_EXEC
#define CODEX_CONT
using CodeXId = std::size_t;
#endif

#endif
