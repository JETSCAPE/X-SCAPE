/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
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

template<class T>
class JetScapeWriterIsrStream : public JetScapeWriterStream<T>
{

 public:

  JetScapeWriterIsrStream<T>() {};
  JetScapeWriterIsrStream<T>(string m_file_name_out):
    JetScapeWriterStream<T>(m_file_name_out){}
  virtual ~JetScapeWriterIsrStream<T>(){};

  //void Init();
  //void Exec();

  void WriteIsr(weak_ptr<PartonShower> ps){
      JetScapeWriterStream<T>::Write(ps);}

  void Write(weak_ptr<PartonShower> ps) {};
  //void Write(weak_ptr<Parton> p) {};
  void Write(weak_ptr<Vertex> v) {
    JetScapeWriterStream<T>::Write(v);
  };
  void Write(weak_ptr<Hadron> h) {};
  //void WriteComment(string s) {};

 private:

};

typedef JetScapeWriterIsrStream<ofstream> JetScapeWriterIsrAscii;
#ifdef USE_GZIP
typedef JetScapeWriterIsrStream<ogzstream> JetScapeWriterIsrAsciiGZ;
#endif

} // end namespace Jetscape

#endif // JETSCAPEWRITERISRSTREAM_H
