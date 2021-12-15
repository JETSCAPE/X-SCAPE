#ifndef CascadeTest_H
#define CascadeTest_H

//#include "Afterburner.h"
#include "BulkMediaBase.h"
#include "FluidCellInfo.h"
#include "BulkMediaInfo.h"
using namespace Jetscape;

class CascadeTest : public BulkMediaBase {
private:
  
  std::vector<std::shared_ptr<Hadron>> hList;

public:

  CascadeTest() {SetId("CascadeTest");};
  virtual ~CascadeTest() {};

  void InitTask();

  virtual void CalculateTime();
  virtual void ExecTime();

  virtual void InitPerEvent();
  virtual void FinishPerEvent();

  virtual any GetHistory(); 

  void GetBulkInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y,
                    Jetscape::real z,
                    std::unique_ptr<BulkMediaInfo> &bulk_info_ptr);
  
};

#endif // CascadeTest_H
