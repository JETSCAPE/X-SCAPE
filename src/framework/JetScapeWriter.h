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

// generic jetscape writer base clase

#ifndef JETSCAPEWRITER_H
#define JETSCAPEWRITER_H

#include <string>
#include "JetScapeModuleBase.h"
#include "PartonShower.h"
#include "JetClass.h"
#include "JetScapeEventHeader.h"

using std::to_string;

namespace Jetscape {

/**
 * @class JetScapeWriter
 * @brief Abstract base class for all JetScape output writers.
 *
 * The `JetScapeWriter` class defines the common interface for all writer
 * modules in the JetScape framework. Writers are responsible for serializing
 * and storing simulation objects, such as partons, hadrons, vertices, and full
 * parton showers.
 *
 * This class also handles writing initialization information and event headers.
 * Specific output formats (ASCII, binary, compressed, etc.) are implemented in
 * derived classes.
 */
class JetScapeWriter : public JetScapeModuleBase {
 public:
  /**
   * @brief Default constructor.
   *
   * Initializes the writer with a default identifier.
   */
  JetScapeWriter() { SetId("JetScape writer"); };

  /**
   * @brief Constructor with output file name.
   * @param m_file_name_out Name of the output file to write.
   */
  JetScapeWriter(string m_file_name_out) { file_name_out = m_file_name_out; }

  /**
   * @brief Destructor.
   */
  virtual ~JetScapeWriter(){};

  /**
   * @brief Set the output file name.
   * @param m_file_name_out Name of the output file.
   */
  void SetOutputFileName(string m_file_name_out) {
    file_name_out = m_file_name_out;
  }

  /**
   * @brief Get the output file name.
   * @return The currently set output file name.
   */
  string GetOutputFileName() { return file_name_out; }

  /**
   * @brief Check if the writer is in a valid state.
   * @return `true` if the writer is active and ready, `false` otherwise.
   */
  virtual bool GetStatus() = 0;

  /**
   * @brief Close the output file/stream.
   */
  virtual void Close(){};

  /**
   * @brief Open the output file/stream.
   */
  virtual void Open(){};

  /**
   * @brief Write initialization XML file.
   */
  virtual void WriteInitFileXML(){};

  /// @name Write functions
  /// Overloaded write methods for various JetScape physics objects.
  /// Derived classes should override the relevant methods to handle output.
  /// @{
  virtual void Write(weak_ptr<Parton> p){};
  virtual void Write(weak_ptr<Jet> j){};
  virtual void Write(weak_ptr<Vertex> v){};
  virtual void Write(weak_ptr<PartonShower> ps){};
  virtual void Write(string s){};
  virtual void WriteComment(string s){};
  virtual void WriteWhiteSpace(string s){};
  virtual void Write(ostream *o){};
  virtual void Write(weak_ptr<Hadron> h){};
  virtual void Write(weak_ptr<Qvector> Qv){};
  /// @}

  /**
   * @brief Write event header to file.
   *
   * Called before event content is written. Can be used to add metadata such as
   * cross section, event weight, number of participants, etc.
   */
  virtual void WriteHeaderToFile(){};

  /**
   * @brief Finalize writing of an event.
   *
   * Called after all modules have written their event data.
   */
  virtual void WriteEvent(){};

  /**
   * @brief Get the event header object.
   * @return Reference to the current `JetScapeEventHeader`.
   */
  virtual JetScapeEventHeader &GetHeader() { return header; };

 protected:
  string file_name_out;        ///< Output file name
  JetScapeEventHeader header;  ///< Event header information
};

}  // end namespace Jetscape

#endif
