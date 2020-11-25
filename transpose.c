#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <omp.h>
#include <pthread.h>

#include <sys/mman.h>

#define ELAPSED_NS(start,stop) \
  (((int64_t)stop.tv_sec-start.tv_sec)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

#define USE_MULTI_THREAD

#define CORRECT_INPUT 10
#define NTHREADS 4


typedef struct
{
  size_t ntime;
  size_t obsnchan;
  size_t nbits;
  size_t ndim;
  size_t npol;
  size_t piperblk;
} db_transpose_t;

/*
int transpose(db_transpose_t * ctx, const void* in, void* out)
{
  // To be used in pointer arithmetic
  const char* inbuf;
  char* outbuf;

  const char* baseinbuf = in;
  char* baseoutbuf      = out;

  // number of packets that span the entire data block, in time
  size_t itime_packets = ctx->piperblk / ctx->ntime;

  // bytes to stride within a packet in input buffer
  // and also amount to copy at a time
  size_t istride = (ctx->npol * ctx->ndim * ctx->nbits * ctx->ntime)/8;

  size_t tstride = ctx->obsnchan * istride;
  size_t ostride = itime_packets * istride;

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
  return 0;
}


int transpose(db_transpose_t * ctx, const void* in, void* out)
{
  // To be used in pointer arithmetic
  const char* inbuf;
  char* outbuf;

  const char* baseinbuf = in;
  char* baseoutbuf      = out;

  // number of packets that span the entire data block, in time
  size_t itime_packets = ctx->piperblk / ctx->ntime;

  // bytes to stride within a packet in input buffer
  // and also amount to copy at a time
  size_t istride = (ctx->npol * ctx->ndim * ctx->nbits * ctx->ntime)/8;

  size_t tstride = ctx->obsnchan * istride;
  size_t ostride = itime_packets * istride;

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
  fprintf(stderr, "transpose...\n");
  return 0;
}
*/

int *transpose_pthread(void* threadid, db_transpose_t * ctx, const void* in, void* out)
{
  fprintf(stderr,"Using pthreads threads\n");
  // To be used in pointer arithmetic
  const char* inbuf;
  char* outbuf;

  const char* baseinbuf = in;
  char* baseoutbuf      = out;

  // number of packets that span the entire data block, in time
  size_t itime_packets = ctx->piperblk / ctx->ntime;

  // bytes to stride within a packet in input buffer
  // and also amount to copy at a time
  //size_t istride = (ctx->npol * ctx->ndim * ctx->nbits * ctx->ntime)/8;
  size_t istride = 32;

  size_t tstride = ctx->obsnchan * istride;
  size_t ostride = itime_packets * istride;

  long tid = (long) threadid;

  // Loop over entire spectrum-packets over all the chans
  for (size_t iptime=tid; iptime < itime_packets; iptime + NTHREADS)
  {
    inbuf  = baseinbuf + iptime*tstride;
    outbuf = baseoutbuf + iptime*istride;
    for (size_t ichan=0; ichan < ctx->obsnchan; ichan+=1)
    {
      memcpy(outbuf, inbuf, istride);
      inbuf += istride;
      outbuf += ostride;
      /*
      for (size_t t=0; t<istride; t++)
        *outbuf++ = *inbuf++;
      outbuf += nstride;
      */
    }
  }
  return 0;
}


int transpose(db_transpose_t * ctx, const void* in, void* out)
{
  // To be used in pointer arithmetic
  const char* inbuf;
  char* outbuf;

  const char* baseinbuf = in;
  char* baseoutbuf      = out;

  // number of packets that span the entire data block, in time
  //size_t itime_packets = ctx->piperblk / ctx->ntime;
  //size_t nchan = ctx->obsnchan;
  size_t itime_packets = 8192*8 / 16;
  size_t nchan = 2048;
  

  // bytes to stride within a packet in input buffer
  // and also amount to copy at a time
  size_t istride = (ctx->npol * ctx->ndim * ctx->nbits * ctx->ntime)/8;

#ifdef USE_MULTI_THREAD
  size_t tstride = ctx->obsnchan * istride;
#endif
  size_t ostride = itime_packets * istride;
  size_t nstride = ostride - istride;

  // Loop over entire spectrum-packets over all the chans
#ifdef USE_MULTI_THREAD
#pragma omp parallel for private (inbuf, outbuf)
#else
  inbuf  = baseinbuf;// + iptime*tstride;
#endif
  for (size_t iptime=0; iptime < itime_packets; iptime++)
  {
#ifdef USE_MULTI_THREAD
    inbuf  = baseinbuf + iptime*tstride;
#endif
    outbuf = baseoutbuf + iptime*istride;
    for (size_t ichan=0; ichan < nchan; ichan+=1)
    {
      memcpy(outbuf, inbuf, istride);
      inbuf += istride;
      outbuf += ostride;
      /*
      for (size_t t=0; t<istride; t++)
        *outbuf++ = *inbuf++;
      outbuf += nstride;
      */
    }
  }
  return 0;
}



/*
int transpose(db_transpose_t * ctx, const void* in, void* out)
{
  // To be used in pointer arithmetic
  const char* inbuf;
  char* outbuf;

  const char* baseinbuf = in;
  char* baseoutbuf      = out;

  // number of packets that span the entire data block, in time
  size_t itime_packets = ctx->piperblk / ctx->ntime;

  // bytes to stride within a packet in input buffer
  // and also amount to copy at a time
  size_t istride = (ctx->npol * ctx->ndim * ctx->nbits * ctx->ntime)/8;

  size_t tstride = ctx->obsnchan * istride;
  size_t ostride = itime_packets * istride;

  // Loop over entire spectrum-packets over all the chans
#pragma omp parallel for private (inbuf, outbuf)
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
  fprintf(stderr, "transpose...\n");
  return 0;
}
*/

int test_buf(char* outbuf, size_t nsamps)
{
  for (size_t i=0; i< nsamps; i++)
  {
    if ( outbuf[i] != CORRECT_INPUT )
    {
      fprintf(stderr, "Test failed at sample: %li, value: %i\n", i, outbuf[i]);
      return EXIT_FAILURE;
    }
  }
  return 0;
}


int main(int argc, char* argv[])
{
  size_t NTIME    = 16; //ntime samples per packet. This is the only constant in the entire pipeline
  size_t OBSNCHAN = 2048; //total number of channels (all the antennas)
  size_t NBITS    = 4; //Number of bits
  size_t NDIM     = 2; //i.e. complex
  size_t NPOL     = 2; //Number of polarisations
  size_t PIPERBLK = 8192*8; // Number of time samples in a block

#ifdef USE_MULTI_THREAD
  omp_set_num_threads(NTHREADS);
#endif

  // create context
  db_transpose_t ctx;
  ctx.ntime = NTIME;
  ctx.obsnchan = OBSNCHAN;
  ctx.nbits = NBITS;
  ctx.ndim = NDIM;
  ctx.npol = NPOL;
  ctx.piperblk = PIPERBLK;
  

  //Create input and output buffers, same size
  size_t buffsize = (NPOL * NDIM * NBITS * PIPERBLK * OBSNCHAN)/8; // Bytes

  fprintf(stderr, "BUFFSIZE: %li\n", buffsize); 

  char* outbuf = malloc(buffsize);
  memset(outbuf, 0, buffsize);
  char* inbuf = malloc(buffsize);
  memset(inbuf, CORRECT_INPUT, buffsize);

  /*
  if (mlock(inbuf, buffsize) == -1)
    fprintf(stderr,"Warning: could not lock memory\n"); 
  if (mlock(outbuf, buffsize) == -1)
    fprintf(stderr,"Warning: could not lock memory\n"); 
  */

  struct timespec start, stop;

  //pthread_t thread[NTHREADS];

  // Time the function
  clock_gettime(CLOCK_MONOTONIC, &start);
  transpose(&ctx, inbuf, outbuf);
  clock_gettime(CLOCK_MONOTONIC, &stop);

  uint64_t elapsed_ns = ELAPSED_NS(start, stop);
  fprintf(stdout, "Fluffed %.1f MBytes in %.2f ms (%.2f Gb/s)\n", buffsize/1.0e6, elapsed_ns / 1.0e6, (8 * (float)buffsize / elapsed_ns));

  clock_gettime(CLOCK_MONOTONIC, &start);
  for (int i=0; i<100; i++)
    transpose(&ctx, inbuf, outbuf);
  clock_gettime(CLOCK_MONOTONIC, &stop);

  elapsed_ns = ELAPSED_NS(start, stop);
  fprintf(stdout, "Fluffed %.1f MBytes in %.2f ms (%.2f Gb/s)\n", buffsize*100/1.0e6, elapsed_ns / 1.0e6, (8 * (float)buffsize*100 / elapsed_ns));

  if (test_buf(outbuf, buffsize) == 0)
    fprintf(stderr, "output buffer passed the test\n");

  free(inbuf);
  free(outbuf);

  return EXIT_SUCCESS;
}

