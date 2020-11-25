#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

#include <sys/mman.h>

#define ELAPSED_NS(start,stop) \
  (((int64_t)stop.tv_sec-start.tv_sec)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

#define USE_MUTLI_THREAD

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

typedef struct
{
  size_t ntime;
  size_t obsnchan;
  size_t nbits;
  size_t ndim;
  size_t npol;
  size_t piperblk;

  size_t itime_packets;
  size_t istride;
  size_t tstride;
  size_t ostride;
} db_transpose_t;


__global__
void transpose(db_transpose_t * ctx, const void* in, void* out)
{
  // To be used in pointer arithmetic
  const char* inbuf;
  char* outbuf;

  const char* baseinbuf = (const char*) in;
  char* baseoutbuf      = (char*) out;

  //size_t itime_packets = ctx->itime_packets;
  size_t istride = ctx->istride;
  size_t tstride = ctx->tstride;
  size_t ostride = ctx->ostride;

  size_t ichan;
  size_t iptime = blockIdx.x;
  size_t nchans = ctx->obsnchan;

  size_t nloops = (ctx->obsnchan + blockDim.x - 1)/blockDim.x;

  for (int igrid=0; igrid < nloops; igrid++)
  {
      ichan = threadIdx.x + blockDim.x * igrid;

      if (ichan < nchans)
      {
          inbuf  = baseinbuf  + iptime*tstride + ichan*istride;
          outbuf = baseoutbuf + iptime*istride + ichan*ostride;
          memcpy(outbuf, inbuf, istride);
      }
  }

  /*
  // Loop over entire spectrum-packets over all the chans
//#pragma omp parallel for private (inbuf, outbuf)
  for (size_t iptime=0; iptime < itime_packets; iptime++)
  {
    inbuf  = baseinbuf + iptime*tstride;
    outbuf = baseoutbuf + iptime*istride;
    for (size_t ichan=0; ichan < ctx->obsnchan; ichan++)
    {
      memcpy(outbuf, inbuf, istride);
      inbuf += istride;
      outbuf += ostride;
    }
  }
  */
}





int test_buf(char* outbuf, size_t nsamps)
{
  for (size_t i=0; i< nsamps; i++)
  {
    if ( outbuf[i] != 10 )
    {
      fprintf(stderr, "Test failed at sample: %li, value: %i\n", i, outbuf[i]);
      return EXIT_FAILURE;
    }
  }
  return 0;
}


int main(int argc, char* argv[])
{
  size_t NTIME    = 16; //ntime samples per packet
  size_t OBSNCHAN = 2048; //total number of channels (all the antennas)
  size_t NBITS    = 4; //Number of bits
  size_t NDIM     = 2; //i.e. complex
  size_t NPOL     = 2; //Number of polarisations
  size_t PIPERBLK = 8192*8; // Number of time samples in a block

  //size_t NTHREADS = 1;
  //omp_set_num_threads(NTHREADS);

  // create context
  db_transpose_t* ctx;
  cudaMallocManaged((void**)&ctx, sizeof(*ctx) );

  ctx->ntime = NTIME;
  ctx->obsnchan = OBSNCHAN;
  ctx->nbits = NBITS;
  ctx->ndim = NDIM;
  ctx->npol = NPOL;
  ctx->piperblk = PIPERBLK;

  // number of packets that span the entire data block, in time
  ctx->itime_packets = ctx->piperblk / ctx->ntime;

  // bytes to stride within a packet in input buffer
  // and also amount to copy at a time
  ctx->istride = (ctx->npol * ctx->ndim * ctx->nbits * ctx->ntime)/8;

  ctx->tstride = ctx->obsnchan * ctx->istride;
  ctx->ostride = ctx->itime_packets * ctx->istride;

  fprintf(stderr, "itime packets: %li, obsnchan: %li\n",
          ctx->itime_packets, ctx->obsnchan);


  //Create input and output buffers, same size
  size_t buffsize = (NPOL * NDIM * NBITS * PIPERBLK * OBSNCHAN)/8; // Bytes

  fprintf(stderr, "BUFFSIZE: %li\n", buffsize); 

  char* outbuf = (char*) malloc(buffsize);
  char* inbuf = (char*) malloc(buffsize);

  for (size_t i=0; i<buffsize; i++)
  {
      outbuf[i] = 0;
      inbuf[i]  = 10;
  }

  // GPU memory
  char *d_inbuf, *d_outbuf;
  if (cudaMalloc((void**) &d_inbuf, buffsize) != cudaSuccess)
  {
      fprintf(stderr, "Couldn't allocate inbuf\n");
      return -1;
  }
  if (cudaMalloc((void**) &d_outbuf, buffsize) != cudaSuccess)
  {
      fprintf(stderr, "Couldn't allocate outbuf\n");
      return -1;
  }

  struct timespec start, stop;

  // Time the function
  clock_gettime(CLOCK_MONOTONIC, &start);
  int nloops = 60;
  for (size_t i = 0; i<nloops; i++)
  {
      if (cudaMemcpy(d_inbuf, inbuf, buffsize, cudaMemcpyHostToDevice) != cudaSuccess)
      {
          fprintf(stderr, "Couldn't move data to device\n");
          return -1;
      }
      //transpose <<<ctx->itime_packets, 512>>> (ctx, d_inbuf, d_outbuf);
      transpose <<<ctx->itime_packets, 1024>>> (ctx, d_inbuf, d_outbuf);
      //cudaMemcpy(d_outbuf, d_inbuf, buffsize, cudaMemcpyDeviceToDevice);
      cudaDeviceSynchronize();

      cudaError_t error = cudaGetLastError();
      if (error != cudaSuccess)
      {
          fprintf(stderr, "Couldn't move data back to host, retval: %s\n", cudaGetErrorString(error));
          return -1;
      }
      cudaMemcpy(outbuf, d_outbuf, buffsize, cudaMemcpyDeviceToHost); 
  }

  clock_gettime(CLOCK_MONOTONIC, &stop);

  uint64_t elapsed_ns = ELAPSED_NS(start, stop);
  fprintf(stdout, "Fluffed %.1f MBytes in %.2f ms (%.2f Gb/s)\n", buffsize*nloops/1.0e6, elapsed_ns / 1.0e6, (8 * (float)buffsize*nloops / elapsed_ns));

  test_buf(outbuf, buffsize);

  free(inbuf);
  free(outbuf);

  cudaFree(d_inbuf);
  cudaFree(d_outbuf);

  return EXIT_SUCCESS;
}

