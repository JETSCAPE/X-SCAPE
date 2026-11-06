/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>

#include <iostream>

using namespace std;

#define MAGENTA "\033[35m"

namespace Jetscape {

HardProcess::HardProcess() {
  VERBOSE(8);
  SetId("HardProcess");
}

HardProcess::~HardProcess() {
  VERBOSE(8);
  hp_list.clear();
  hd_list.clear();
  ps_list.clear();

  disconnect_all();
}

void HardProcess::Init() {
  JetScapeModuleBase::InitTask();

  JSINFO << "Initialize HardProcess : " << GetId() << " ...";

  VERBOSE(8);

  ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
  if (!ini) {
    // If not vacuum case, give warning to add initial state module
    bool in_vac = GetXMLElementInt({"Eloss", "Matter", "in_vac"});
    bool in_brick = GetXMLElementInt({"Eloss", "Matter", "brick_med"});
    if (!in_vac and !in_brick) {
      JSWARN << "No initial state module! Please check whether you intend to "
                "add an initial state module.";
      exit(-1);
    }
  }
  string status = GetXMLElementText({"PartonPrinter", "Status"});
  if (status != "off") {
    printer = GetXMLElementText({"PartonPrinter", "FileName"});
    JSINFO << BOLDYELLOW << "Extra parton info goes to " << printer;
  }

  InitTask();
  InitTasks();
}

void HardProcess::ExecuteTask() {
  JSINFO << "Run Hard Process : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();

  // Subtasks assumed here to be ISR related ...
}

void HardProcess::ClearTask() {
  JSDEBUG << "Clear Hard Process : " << GetId() << " ...";

  hp_list.clear();
  hd_list.clear();
  ps_list.clear();
  Remnants.clear();
  color_tag_ = 0;
  extra_color_tag_.clear();
  block_color_tag_.clear();
  block_color_base_ = 0;

  VERBOSE(8) << hp_list.size();
}

unsigned int HardProcess::NextColorTag(int scatter, int block) {
  if (block >= 0) {
    if (block >= static_cast<int>(block_color_tag_.size()))
      throw std::runtime_error("HardProcess::NextColorTag: colour tag block " +
                               std::to_string(block) +
                               " was not reserved for this event");
    unsigned int &tag = block_color_tag_[block];
    if (tag + 1 >= block_color_base_ + (block + 1u) * color_block_size)
      throw std::runtime_error(
          "HardProcess::NextColorTag: a shower used up its " +
          std::to_string(color_block_size) + " colour tags");
    return ++tag;
  }
  if (scatter > 0) {
    unsigned int &tag = extra_color_tag_[scatter];
    const unsigned int floor = 100000u * scatter + 50000u;
    if (tag < floor)
      tag = floor;
    return ++tag;
  }
  unsigned int n = 0;  // the triggering scattering's hard partons
  for (int i = 0; i < static_cast<int>(hp_list.size()); ++i)
    n += GetHardPartonScatter(i) == 0;
  const unsigned int floor = max_colorPerShower * (n + 1);
  if (color_tag_ < floor)
    color_tag_ = floor;
  return ++color_tag_;
}

void HardProcess::ReserveColorTagBlocks(int n) {
  // Above the counters of NextColorTag (and the tags they start above): the
  // triggering scattering's, from max_colorPerShower x (hard partons + 1), and
  // the additional scatterings', below 100000 x (scatter + 1); and above every
  // tag the hard partons carry and every tag handed out so far.
  unsigned int top = max_colorPerShower * (hp_list.size() + 1) + color_tag_;
  for (int i = 0; i < static_cast<int>(hp_list.size()); ++i) {
    top = std::max(top, 100000u * (GetHardPartonScatter(i) + 1u));
    top = std::max({top, hp_list[i]->color(), hp_list[i]->anti_color()});
  }
  for (const auto &kv : extra_color_tag_)
    top = std::max(top, kv.second);
  const unsigned long long base =
      (top / color_block_size + 1ull) * color_block_size;
  const unsigned long long end =
      base + static_cast<unsigned long long>(n) * color_block_size;
  if (end > static_cast<unsigned long long>(std::numeric_limits<int>::max()))
    throw std::runtime_error(
        "HardProcess::ReserveColorTagBlocks: " + std::to_string(n) +
        " blocks of colour tags do not fit below 2^31");
  block_color_base_ = static_cast<unsigned int>(base);
  block_color_tag_.resize(n);
  for (int b = 0; b < n; ++b)
    block_color_tag_[b] = block_color_base_ + b * color_block_size;
}

void HardProcess::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);

  auto f = w.lock();
  if (f) {
    VERBOSE(8) << f->GetOutputFileName();

    // Weight, xsec, etc

    // // Can explicitly write our own header information, though the writer
    // should handle this. std::ostringstream oss; oss.str(""); oss << GetId()
    // << " sigmaGen  = " << GetSigmaGen(); f->WriteComment ( oss.str() );
    // oss.str(""); oss << GetId() << " sigmaErr  = " << GetSigmaErr();
    // f->WriteComment ( oss.str() );
    // oss.str(""); oss << GetId() << " weight  = " << GetEventWeight();
    // f->WriteComment ( oss.str() );

    // Hard partons
    f->WriteComment("HardProcess Parton List: " + GetId());
    for (auto hp : hp_list)
      f->Write(hp);

    // Commented out for now, decide/fix how to safe store the ISR shower ...
    // JetScapeTask::WriteTasks(w);
  }
}

void HardProcess::CollectHeader(weak_ptr<JetScapeWriter> w) {
  auto f = w.lock();
  if (f) {
    auto &header = f->GetHeader();
    header.SetSigmaGen(GetSigmaGen());
    header.SetSigmaErr(GetSigmaErr());
    header.SetPtHat(GetPtHat());
    header.SetEventWeight(GetEventWeight());

    if (hp_list.size() > 0) {
      header.SetVertexX(hp_list[0]->x_in().x());
      header.SetVertexY(hp_list[0]->x_in().y());
      header.SetVertexZ(hp_list[0]->x_in().z());
    }
  }
}

}  // end namespace Jetscape
