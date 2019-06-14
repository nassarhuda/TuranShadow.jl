#include "stdafx.h"

#include "n2v.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

void LoadEdges(int64* edges, int64 nedges, bool& Directed, PWNet& InNet)
{
  for (int64 ei = 0; ei < nedges; ++ei) {
    int64 SrcNId = edges[2*ei];
    int64 DstNId = edges[2*ei+1];
    //printf("%i %i\n", (int)SrcNId, (int)DstNId);

    double Weight = 1.0;
    if (!InNet->IsNode(SrcNId)){ InNet->AddNode(SrcNId); }
    if (!InNet->IsNode(DstNId)){ InNet->AddNode(DstNId); }
    InNet->AddEdge(SrcNId,DstNId,Weight);
    if (!Directed){ InNet->AddEdge(DstNId,SrcNId,Weight); }
  }
}

void StoreOutput(int64 nnodes, int64* outids, double* embed,
    TIntFltVH& EmbeddingsHV) {
  size_t oi = 0;
  size_t ei = 0;
  for (int i = EmbeddingsHV.FFirstKeyId(); EmbeddingsHV.FNextKeyId(i);) {
    outids[oi++] = EmbeddingsHV.GetKey(i);
    for (int64 j = 0; j < EmbeddingsHV[i].Len(); j++) {
      embed[ei++] = EmbeddingsHV[i][j];
    }
  }
}

extern "C" {
int node2vec(int64 *edges, int64 nedges, int64 nnodes, // input
             int64 *outids, double *outemb,            // output
             int Dimensions, int WalkLen, int NumWalks,
             int WinSize, int Iter,
             double ParamP, double ParamQ) // parameters
{
  //int Dimensions=d, WalkLen=80, NumWalks=10, WinSize=10, Iter=1;
  //double ParamP=1.0, ParamQ=1.0;
  bool Directed=false, Weighted=false, Verbose=false;

  PWNet InNet = PWNet::New();
  TIntFltVH EmbeddingsHV;

  //printf("%i\n", (int)nedges);

  LoadEdges(edges, nedges, Directed, InNet);

  node2vec(InNet, ParamP, ParamQ, Dimensions, WalkLen, NumWalks, WinSize, Iter,
   Verbose, EmbeddingsHV);

  StoreOutput(nnodes, outids, outemb, EmbeddingsHV);

  return 1;
}

} // extern "C"
