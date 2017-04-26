#include "algebra/graph/bgmgraph.h"

namespace pascinference {
namespace algebra {

extern void cuda_BGMGraph_destroy(int *neighbor_nmbs_gpu, int **neighbor_ids_gpu, int **neighbor_ids_cpugpu, int n);
extern void cuda_BGMGraph_process(int *neighbor_nmbs_gpu, int **neighbor_ids_gpu, int **neighbor_ids_cpugpu, int n);

}
} /* end of namespace */
