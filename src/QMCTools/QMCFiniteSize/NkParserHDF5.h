#ifndef NK_PARSER_HDF5_H
#define NK_PARSER_HDF5_H

#include "QMCTools/QMCFiniteSize/NkParserBase.h"
#include "io/hdf/hdf_archive.h"

namespace qmcplusplus
{
class NkParserHDF5 : public NkParserBase
{
public:
  void parse(const string& fname) override;

private:
  hdf_archive statfile;
};
} // namespace qmcplusplus
#endif
