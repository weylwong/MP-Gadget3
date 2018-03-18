#ifndef GENIC_PROTO_H
#define GENIC_PROTO_H
#include <bigfile.h>
#include <stdint.h>

/* shift is the shift that will be applied when saving, in units of a grid cell.
 * It is a 3-vector, each element corresponding to an axis.*/
void   displacement_fields(int Ngrid, int Type, const double* shift);

uint64_t id_offset_from_index(const int i, const int Ngrid);
void   free_ffts(void);

void saveheader(BigFile * bf, int64_t TotNumPart, int64_t TotNuPart, double nufrac);
void  write_particle_data(const int Type, BigFile * bf,  const uint64_t FirstID, const double *shift, const int Ngrid);

void  read_parameterfile(char *fname);
#endif
