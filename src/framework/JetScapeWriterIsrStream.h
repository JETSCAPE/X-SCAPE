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

#ifndef JETSCAPEWRITERISRSTREAM_H
#define JETSCAPEWRITERISRSTREAM_H

#include "JetScapeWriterStream.h"

namespace Jetscape {

/**
 * @file JetScapeWriterIsrStream.h
 * @brief Specialised writer stream for ISR-related objects.
 *
 * This header defines a small template writer that forwards selected
 * writes to the base `JetScapeWriterStream<T>` implementation. The
 * template parameter `T` is the underlying output stream type (for
 * example `std::ofstream` or a gzipped stream type when `USE_GZIP` is
 * defined).
 */
template <class T>
class JetScapeWriterIsrStream : public JetScapeWriterStream<T> {
 public:
  /**
   * @brief Default constructor.
   *
   * Constructs an instance without opening a file. Use the other
   * constructor or the base class API to open a file for writing.
   */
  JetScapeWriterIsrStream<T>(){};

  /**
   * @brief Construct and open an output file.
   * @param m_file_name_out The path to the output file to open.
   *
   * Forwards the filename to the base `JetScapeWriterStream<T>` which
   * handles opening the stream.
   */
  JetScapeWriterIsrStream<T>(string m_file_name_out)
      : JetScapeWriterStream<T>(m_file_name_out) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~JetScapeWriterIsrStream<T>(){};

  // void InitTask();
  // void ExecuteTask();

  /**
   * @brief Write an ISR `PartonShower` to the stream.
   * @param ps Weak pointer to the `PartonShower` object to write.
   *
   * This helper calls the base class `Write` implementation to serialize
   * the provided `PartonShower` into the underlying stream.
   */
  void WriteIsr(weak_ptr<PartonShower> ps) {
    JetScapeWriterStream<T>::Write(ps);
  }

  /**
   * @brief No-op override for `PartonShower` writes.
   * @param ps Weak pointer to a `PartonShower`.
   *
   * Present as an override to suppress or change default behaviour in
   * specific writer specialisations. Intentionally left empty.
   */
  void Write(weak_ptr<PartonShower> ps){};
  // void Write(weak_ptr<Parton> p) {};
  /**
   * @brief Write a `Vertex` to the stream.
   * @param v Weak pointer to the `Vertex` object to write.
   *
   * Forwards to the base class `Write` to perform serialization.
   */
  void Write(weak_ptr<Vertex> v) { JetScapeWriterStream<T>::Write(v); };

  /**
   * @brief No-op override for `Hadron` writes.
   * @param h Weak pointer to a `Hadron`.
   *
   * Included for API completeness; this writer does not output hadrons.
   */
  void Write(weak_ptr<Hadron> h){};
  // void WriteComment(string s) {};

 private:
};

typedef JetScapeWriterIsrStream<ofstream> JetScapeWriterIsrAscii;
#ifdef USE_GZIP
typedef JetScapeWriterIsrStream<ogzstream> JetScapeWriterIsrAsciiGZ;
#endif

}  // end namespace Jetscape

#endif  // JETSCAPEWRITERISRSTREAM_H
