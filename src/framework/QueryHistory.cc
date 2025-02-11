#include "QueryHistory.h"
#include "JetScapeLogger.h"
#include <stdlib.h>
#include <algorithm>
#include <memory>

using namespace std;

namespace Jetscape {

QueryHistory *QueryHistory::m_pInstance = NULL;

QueryHistory *QueryHistory::Instance() {
  if (!m_pInstance) {
    JSINFO << "Created QueryHistory Instance";
    m_pInstance = new QueryHistory();
  }

  return m_pInstance;
}

any QueryHistory::GetHistoryFromModule(string mName)
{
  //JP: TO be implement and or only use the FromModules ...
  auto it = taskMap.find(mName);

  if (it == taskMap.end()) {
    JSWARN << "Module name " << mName << " not found in taskMap";
    return 0;
  }

  auto module_ptr = it->second.lock();
  if (!module_ptr) {
    JSWARN << "Weak pointer to module " << mName << " is expired";
    return 0;
  }

  auto history_ptr = std::dynamic_pointer_cast<JetScapeModuleBase>(module_ptr);
  if (!history_ptr) {
    JSWARN << "Unable to cast to JetScapeModuleBase for module " << mName;
    return 0;
  }

  return history_ptr->GetHistory();
}

vector<any> QueryHistory::GetHistoryFromModules(string mName)
{
  int num = taskMap.count(mName);
  //DEBUG:
  //cout<<"--> "<<num<<endl;
  vector<any> mHistories;
  auto it = taskMap.equal_range(mName);

  for (auto itr = it.first; itr != it.second; ++itr)
  {
    //JP: maybe change JetScapeTask -> JetScapeModuleBase in header to avoid this dynamic casting etc to be followed up ...
    if (std::dynamic_pointer_cast<JetScapeModuleBase>(itr->second.lock()))
      mHistories.push_back(std::dynamic_pointer_cast<JetScapeModuleBase>(itr->second.lock())->GetHistory());
  }

  return mHistories;
}

void QueryHistory::UpdateTaskMap()
{

  VERBOSE(2) << "QueryHistory::UpdateTaskMap()";

  //JP: Think about smarter/more efficient way rather than clear map and iterate through all tasks again ...
  taskMap.clear();
  auto mt = main_task.lock();

  //Quick and dirty to see all tasks ... make recursive if needed
  if (mt) {

    for (auto it : mt->GetTaskList())
    {

      //JSINFO << it->GetId();
      taskMap.emplace(it->GetId(), it);

      for (auto it2 : it->GetTaskList())
      {
        //JSINFO  << " " << it2->GetId() ;
        taskMap.emplace(it2->GetId(), it2);
      }
    }
  }
}

void QueryHistory::PrintTaskMap()
{
  JSINFO << "QueryHistory::PrintTaskMap()";

  for (auto& x : taskMap) {
    JSINFO << " " << x.first << ":\t " << x.second.lock().get() << "\t active = " << x.second.lock()->GetActive() << "\t multiThread = " << x.second.lock()->GetMultiThread()<<"\t Task number = "<<x.second.lock()->GetMyTaskNumber();
    if (std::dynamic_pointer_cast<JetScapeModuleBase>(x.second.lock())) JSINFO<<" \t \t \t IsTimeStepped = "<<std::dynamic_pointer_cast<JetScapeModuleBase>(x.second.lock())->IsTimeStepped();
  }
}

void QueryHistory::PrintTasks()
{
  //Quick and dirty to see all tasks ... make recursive ...

  JSINFO << "QueryHistory::PrintTasks()";

  auto mt = main_task.lock();

  //Quick and dirty to see all tasks ... make recursive ...
  if (mt) {
    for (auto it : mt->GetTaskList()) {
      JSINFO << it->GetId();
      for (auto it2 : it->GetTaskList()) {
        JSINFO  << " " << it2->GetId() ;
        for (auto it3 : it2->GetTaskList())
          JSINFO  << "  " << it3->GetId() ;
      }
    }
  }
}

}
