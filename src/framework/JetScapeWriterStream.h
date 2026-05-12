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

// jetscape writer ascii class

#ifndef JETSCAPEWRITERSTREAM_H
#define JETSCAPEWRITERSTREAM_H

#include <fstream>
#include <string>

#ifdef USE_GZIP
#include "gzstream.h"
#endif

#include "JetScapeWriter.h"

using std::ofstream;

namespace Jetscape {

/**
 * @class JetScapeWriterStream
 * @brief Writer module for exporting JETSCAPE events to an ASCII (or gzipped)
 * stream.
 *
 * The JetScapeWriterStream class provides functionality for writing
 * simulation outputs such as partons, vertices, hadrons, and entire
 * parton showers into human-readable ASCII files.
 *
 * The implementation supports both plain-text (`ofstream`) and
 * compressed (`ogzstream`) backends.
 *
 * @tparam T Output stream type (`ofstream` or `ogzstream`).
 */
template <class T>
class JetScapeWriterStream : public JetScapeWriter {
 public:
  /**
   * @brief Default constructor.
   */
  JetScapeWriterStream<T>(){};

  /**
   * @brief Constructor with output file name.
   * @param m_file_name_out Name of the output file.
   */
  JetScapeWriterStream<T>(string m_file_name_out);

  /**
   * @brief Destructor.
   */
  virtual ~JetScapeWriterStream<T>();

  /** @brief Initialize the writer module.
   *
   * Opens the output file.
   */
  void InitTask();

  /**
   * @brief Execute the writer for an event.
   *
   * Normally invoked in the event loop. By default,
   * this does nothing, since writing is triggered
   * by module calls to `Write()`.
   */
  void ExecuteTask();

  /**
   * @brief Get the output file status.
   * @return True if the output file stream is in a good state.
   */
  bool GetStatus() { return output_file.good(); }

  /**
   * @brief Close the output file stream.
   */
  void Close() { output_file.close(); }

  /**
   * @brief Write the main XML configuration file into the output stream.
   */
  void WriteInitFileXMLMain();

  /**
   * @brief Write the user XML configuration file into the output stream.
   */
  void WriteInitFileXMLUser();

  /**
   * @brief Write a parton shower to the output file.
   * @param ps Weak pointer to a PartonShower object.
   *
   * Outputs vertices and partons belonging to the shower in a structured
   * format.
   */
  void Write(weak_ptr<PartonShower> ps);

  /**
   * @brief Write a single parton to the output file.
   * @param p Weak pointer to a Parton object.
   */
  void Write(weak_ptr<Parton> p);

  /**
   * @brief Write a vertex to the output file.
   * @param v Weak pointer to a Vertex object.
   */
  void Write(weak_ptr<Vertex> v);

  /**
   * @brief Write a hadron to the output file.
   * @param h Weak pointer to a Hadron object.
   */
  void Write(weak_ptr<Hadron> h);
  // void Write(weak_ptr<Qvector> Qv);

  /**
   * @brief Write event header metadata to the output file.
   *
   * Includes information such as:
   * - Cross section and error.
   * - Event weight.
   * - Number of participants (`Npart`).
   * - Number of binary collisions (`Ncoll`).
   * - Total entropy.
   * - Event plane angle.
   */
  void WriteHeaderToFile();

  /**
   * @brief Write a raw string to the output file.
   * @param s The string to write.
   */
  void Write(string s) { output_file << s << endl; }

  /**
   * @brief Write a comment line to the output file.
   * @param s The string to write as a comment (prefixed with '#').
   */
  void WriteComment(string s) { output_file << "# " << s << endl; }

  /**
   * @brief Write whitespace-separated tokens to the output file.
   * @param s The token string.
   */
  void WriteWhiteSpace(string s) { output_file << s << " "; }

  /**
   * @brief Write an event marker.
   *
   * By default, does nothing. Event writing is handled by the framework.
   */
  void WriteEvent();

 protected:
  T output_file;  //!< Output file stream (ASCII or gzip-compressed).
  int hadronCounter = 0;
  // int m_precision; //!< Output precision

  // Allows the registration of the module so that it is available to be used by
  // the Jetscape framework.

  /**
   * @brief Register the ASCII writer with the JETSCAPE framework.
   */
  static RegisterJetScapeModule<JetScapeWriterStream<ofstream>> reg;

  /**
   * @brief Register the gzip writer with the JETSCAPE framework (if enabled).
   */
  static RegisterJetScapeModule<JetScapeWriterStream<ogzstream>> regGZ;
};

/**
 * @typedef JetScapeWriterAscii
 * @brief Type alias for ASCII writer using `std::ofstream`.
 */
typedef JetScapeWriterStream<ofstream> JetScapeWriterAscii;

#ifdef USE_GZIP
/**
 * @typedef JetScapeWriterAsciiGZ
 * @brief Type alias for gzip-compressed writer using `ogzstream`.
 */
typedef JetScapeWriterStream<ogzstream> JetScapeWriterAsciiGZ;
#endif

}  // end namespace Jetscape

#endif  // JETSCAPEWRITERSTREAM_H
