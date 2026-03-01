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

// jetscape writer ascii class, filter

#ifndef JETSCAPEWRITERSTREAMFILTER_H
#define JETSCAPEWRITERSTREAMFILTER_H

#include "JetScapeWriterStream.h"
#include "sigslot.h"

/// @name Filter Flags
/// Bitmask flags that define which JetScape objects should be written to the
/// output. These can be combined with bitwise OR to select multiple object
/// types.
/// @{
#define JETSCAPEWRITER_PARTONSHOWER 1
#define JETSCAPEWRITER_PARTON 2
#define JETSCAPEWRITER_VERTEX 4
#define JETSCAPEWRITER_HADRON 8
/// @}

namespace Jetscape {

/**
 * @class JetScapeWriterStreamFilter
 * @brief Filtered stream writer for JetScape output.
 *
 * This templated writer class extends `JetScapeWriterStream` by introducing
 * an output filter based on object type. Only objects selected by the
 * `displayFilter` bitmask will be written to the output file.
 *
 * @tparam T The output stream type (e.g., `ofstream` for plain ASCII,
 *           or `ogzstream` if gzip compression is enabled).
 *
 * Example usage:
 * @code
 *   // Write only hadrons and partons
 *   auto writer = make_shared<JetScapeWriterAsciiFilter>("output.dat",
 *                        JETSCAPEWRITER_HADRON | JETSCAPEWRITER_PARTON);
 * @endcode
 */
template <class T>
class JetScapeWriterStreamFilter : public JetScapeWriterStream<T> {
 public:
  /**
   * @brief Default constructor.
   */
  JetScapeWriterStreamFilter<T>(){};

  /**
   * @brief Constructor with output file name and filter mask.
   *
   * @param m_file_name_out Output file name.
   * @param filter Bitmask specifying which objects to write
   *               (see filter flags above).
   */
  JetScapeWriterStreamFilter<T>(string m_file_name_out, unsigned char filter)
      : JetScapeWriterStream<T>(m_file_name_out) {
    displayFilter = filter;
  }

  /**
   * @brief Destructor.
   */
  virtual ~JetScapeWriterStreamFilter<T>(){};

  /**
   * @brief Signal for retrieving the list of hadrons.
   *
   * This allows external modules to connect and obtain the current hadron list.
   */
  sigslot::signal1<vector<shared_ptr<Hadron>>&> GetHadronList;

  // void InitTask();
  // void ExecuteTask();
  //
  // bool GetStatus() {return output_file.good();}
  // void Close() {output_file.close();}
  //
  // void WriteInitFileXML();

  /**
   * @brief Write a parton shower if the filter allows it.
   * @param ps Weak pointer to the parton shower.
   */
  void Write(weak_ptr<PartonShower> ps) {
    if (displayFilter & JETSCAPEWRITER_PARTONSHOWER) {
      JetScapeWriterStream<T>::Write(ps);
    }
  }

  /**
   * @brief Write a parton if the filter allows it.
   * @param p Weak pointer to the parton.
   */
  void Write(weak_ptr<Parton> p) {
    if (displayFilter & JETSCAPEWRITER_PARTON) {
      JetScapeWriterStream<T>::Write(p);
    }
  }

  /**
   * @brief Write a vertex if the filter allows it.
   * @param v Weak pointer to the vertex.
   */
  void Write(weak_ptr<Vertex> v) {
    if (displayFilter & JETSCAPEWRITER_VERTEX) {
      JetScapeWriterStream<T>::Write(v);
    }
  }

  /**
   * @brief Write a hadron if the filter allows it.
   * @param h Weak pointer to the hadron.
   */
  void Write(weak_ptr<Hadron> h) {
    if (displayFilter & JETSCAPEWRITER_HADRON) {
      JetScapeWriterStream<T>::Write(h);
    }
  }

  // void WriteHeaderToFile();
  //
  // void Write(string s) {output_file<<s<<endl;}
  // void WriteComment(string s) {output_file<<"# "<<s<<endl;}
  // void WriteWhiteSpace(string s) {output_file<<s<<" ";}
  // void WriteEvent();

 protected:
  unsigned char displayFilter;  ///< Bitmask specifying which objects to output.
};

/**
 * @typedef JetScapeWriterAsciiFilter
 * @brief Convenience typedef for ASCII output writer with filtering.
 */
typedef JetScapeWriterStreamFilter<ofstream> JetScapeWriterAsciiFilter;

#ifdef USE_GZIP
/**
 * @typedef JetScapeWriterAsciiGZFilter
 * @brief Convenience typedef for gzip-compressed ASCII output writer with
 * filtering.
 */
typedef JetScapeWriterStreamFilter<ogzstream> JetScapeWriterAsciiGZFilter;
#endif

}  // end namespace Jetscape

#endif  // JETSCAPEWRITERSTREAM_H
