#MESSAGE(${GPU_SOURCE_DIR})
#include_directories(${COMMON_INCLUDES} ${COMMON_CUDA_INCLUDES} )

include_directories(${COMMON_INCLUDES} ${GOOGLETEST_DIR}/include/)

MESSAGE(${GOOGLETEST_DIR}/include/)

CUDA_ADD_LIBRARY(SceGPU
	Solver.h
	Solver.cu
	SceNodes.h
	SceNodes.cu
	SceECM.h
    SceECM.cu
	SceCells.h
	SceCells.cu
	SimulationDomainGPU.h
	SimulationDomainGPU.cu
)

