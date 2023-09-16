#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdarg.h>
#include <cstdint>

#include "runhfp.h"
#include "jmesh.h"
#include "fittingPrimitives.h"


FILE _iob[] = { *stdin, *stdout, *stderr };

extern "C" FILE * __cdecl __iob_func(void)
{
    return _iob;
}

List* fit_hfp(Triangulation* t) {
    unsigned char wtf = 0;
    wtf |= HFP_FIT_PLANES;
    wtf |= HFP_FIT_SPHERES;
    wtf |= HFP_FIT_CYLINDERS;
    HFP_Action hfp(t);
    List* coll = hfp.fit(wtf);
    Triangle* f;
    Node* n;
    int i = 0; FOREACHVTTRIANGLE((&(t->T)), f, n) f->info = (void*)i++;

    return coll;
}

std::vector<float_t> process_hfp(List* coll, int numclusters, const int numtris) {
    Node* n;
    const int tris = numtris;
    static Node* last_node;
    Triangle* t;
    List* tlist, * collapses = coll;

    int index;
    int32_t* idx = new int32_t[numtris];
    last_node = collapses->head();

    if (numclusters == 1)
    {
        index = reinterpret_cast<intptr_t>(last_node->next()->data);
        for (unsigned int i = 0; i < numtris; i++) idx[i] = index;
    }
    else
    {
        for (unsigned int i = 0; i < numtris; i++) idx[i] = 0;
        for (unsigned int i = 0; i < numclusters-1; i++)
        {
            tlist = (List*)last_node->data;
            index = reinterpret_cast<intptr_t>(last_node->next()->next()->data);
            FOREACHVTTRIANGLE(tlist, t, n) idx[reinterpret_cast<intptr_t>(t->info)] = index;
            last_node = (last_node->next() == NULL) ? (NULL) : (last_node->next()->next()->next());
        }
    }

    std::vector<float_t> out_arr(numtris);

    for (int i = 0; i < numtris; i++) {
        out_arr[i] = idx[i];
    }

    return out_arr;
}

namespace py = pybind11;

std::vector<float_t> wrapper_hfp(py::array_t<float> vs, py::array_t<int> fs, int numclusters) {
    std::vector<float> vector = std::vector<float>(vs.data(), vs.data() + vs.size());
    py::buffer_info buf_vs = vs.request(), buf_fs = fs.request();

    float *ptr_vs = static_cast<float *>(buf_vs.ptr);
    int *ptr_fs = static_cast<int *>(buf_fs.ptr);

    List* coll;
    std::vector<float_t> idx;
    Triangulation* nt = new Triangulation;

    // load numpy arrays into Triangulation
    int ioerr = nt->loadNUMPY(ptr_vs, ptr_fs, buf_vs.shape[0], buf_fs.shape[0]);

    if (nt->shells() == 0)
    {
        JMesh::warning("No triangles loaded!\n");
        delete nt;
    }
    else if (nt->shells() > 1)
    {
        JMesh::warning("Only single component meshes are supported!\n"
            "Removing all the smallest components.");
        nt->removeSmallestComponents();
        coll = fit_hfp(nt);
    }
    else coll = fit_hfp(nt);

    // link triangles to certain indices -> try it works with the process_hfp function
    // when also returning this Triangles and Vectors from nt.
    Triangle* f;
    Node* n;
    int i=0; FOREACHVTTRIANGLE((&(nt->T)), f, n) f->info = (void *)i++;

    const int numtris = nt->T.numels();
    idx = process_hfp(coll, numclusters, numtris);

    return idx;

}


PYBIND11_MODULE(hfp, m) {
    m.doc() = R"pbdoc(
        Hierarchical Fitting Primitives
        -----------------------
        .. currentmodule:: hfp
        .. autosummary::
           :toctree: _generate
    )pbdoc";
    m.def("run_hfp", &wrapper_hfp, R"pbdoc(
        Runs HFP
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}


