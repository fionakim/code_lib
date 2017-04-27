"""Tests alternate OBOReader."""

import sys
import timeit
import datetime

# Test local version of goatools
sys.path.insert(0, '..')
from goatools.obo_parser import GODag

#################################################################
# Sub-routines to tests
#################################################################
def test_write_hier_all(name, go_id, dag, out):
    """Prints the entire mini GO hierarchy, with counts of children."""
    out.write('\nTEST {NAME} {GO}: Print all hierarchies:\n'.format(
        NAME=name, GO=go_id))
    dag.write_hier(go_id, out, num_child=True)
    out.write("GOTerm: {}\n".format(dag[go_id].__repr__()))


def test_write_hier_norep(name, go_id, dag, out):
    """Shortens hierarchy report by only printing branches once.

         Prints the 'entire hierarchy' of GO:0000005 the 1st time seen:

           ---     1 GO:0000005    L-02    D-02
           ----     0 GO:0000010   L-03    D-04

         Prints just GO:0000005 (ommit child GO:10) the 2nd time seen:

           ===     1 GO:0000005    L-02    D-02

         '=' is used in hierarchy mark to indicate that the pathes
             below the marked term have already been printed.
    """
    out.write('\nTEST {NAME} {GO}: Print branches just once:\n'.format(
        NAME=name, GO=go_id))
    dag.write_hier(go_id, out, num_child=True, short_prt=True)
    out.write("GOTerm: {}\n".format(dag[go_id]))


def test_write_hier_lim(dag, out):
    """Limits hierarchy list to GO Terms specified by user."""
    go_omit = ['GO:0000005', 'GO:0000010']
    go_ids = [go_id for go_id in dag if go_id not in go_omit]
    out.write('\nTEST {NAME} OMIT: 05->10:\n'.format(NAME="MINI"))
    dag.write_hier("GO:0000001", out, include_only=go_ids)
      #go_marks=[oGO.id for oGO in oGOs_in_cluster])

def test_write_hier_mrk(dag, out):
    """Print all paths, but mark GO Terms of interest. """
    mark_lst = ['GO:0000001', 'GO:0000003', 'GO:0000006', 'GO:0000008', 'GO:0000009']
    out.write('\nTEST {NAME} MARK: 01->03->06->08->09:\n'.format(NAME="MINI"))
    dag.write_hier("GO:0000001", out, go_marks=mark_lst)
      #go_marks=[oGO.id for oGO in oGOs_in_cluster])

def _load_dag(dag_fin, opt_fields=None, out=None):
    """Run numerous tests for various reports."""
    tic = timeit.default_timer()
    dag = GODag(dag_fin, opt_fields)
    toc = timeit.default_timer()
    msg = "Elapsed HMS for OBO DAG load: {}\n\n".format(str(datetime.timedelta(seconds=(toc-tic))))
    if out is not None:
        out.write(msg)
    else:
        sys.stdout.write(msg)
    return dag

#################################################################
# Tests
#################################################################
def test_all(fout=None):
    """Run all tests in this module."""
    out = sys.stdout if fout is None else open(fout, 'w')
    test_full(out)
    test_full(out, opt_fields='def')
    test_full(out, opt_fields=['def'])
    test_full(out, opt_fields=['def', 'xref'])
    test_mini(out)
    # User can specify either 'def' or 'defn' to indicate 'definition'
    test_mini_contents(out, opt_fields='def')
    test_mini_contents(out, opt_fields='defn')
    test_mini_contents(out, opt_fields='xref')
    test_mini_contents(out, opt_fields=['def', 'xref'])
    test_mini_contents(out, opt_fields=['defn', 'xref'])
    test_full_contents(out)
    if fout is not None:
        out.close()
        sys.stdout.write("  WROTE: {}\n".format(fout))

def test_full_contents(out):
    """Ensure that obo file can be read with all optional attributes."""
    all_fields = [
        'comment', 'consider', 'defn', 'is_class_level', 'is_metadata_tag', 'is_transitive',
        'relationship', 'replaced_by', 'subset', 'synonym', 'transitive_over', 'xref']
    dag_all = _load_dag("./go-basic.obo", all_fields, out)
    dag = _load_dag("./go-basic.obo", None, out)
    if len(dag_all) != len(dag):
        raise Exception("FAILED: test_full_contents")

def test_full(out, opt_fields=None):
    """Use OBOReader in default operation."""
    dag_fin = "./go-basic.obo"
    dag = _load_dag(dag_fin, opt_fields, out)
    test_write_hier_all("FULL", "GO:0000009", dag, out)
    test_write_hier_norep("FULL", "GO:0000010", dag, out)

def test_mini(out, opt_fields=None):
    """Run numerous tests for various reports."""
    dag = _load_dag("./data/mini_obo.obo", opt_fields, out)
    test_write_hier_lim(dag, out)
    test_write_hier_mrk(dag, out)
    test_write_hier_all("MINI", "GO:0000009", dag, out)
    test_write_hier_norep("MINI", "GO:0000010", dag, out)

def test_mini_contents(out, opt_fields):
    """Test that optional terms were loaded."""
    dag = _load_dag("./data/mini_obo.obo", opt_fields, out)
    go_1 = dag["GO:0000001"]
    go_2 = dag["GO:0000002"]
    # CHECK 'xref'
    if 'xref' in opt_fields:
        if len(go_1.xref) != 3:
            raise Exception("XREF TEST FAILED")
        if len(go_2.xref) != 1:
            raise Exception("XREF TEST FAILED")
    else:
        if hasattr(go_1, 'xref'):
            raise Exception("XREF TEST FAILED")
        if hasattr(go_2, 'xref'):
            raise Exception("XREF TEST FAILED")
    # CHECK 'def'
    if 'def' in opt_fields or 'defn' in opt_fields:
        # 'def' is a keyword in obo, but a reserved word in Python,
        # so access 'def' items in an obo using 'defn' in GODag.
        if not isinstance(go_1.defn, str):
            raise Exception("def TEST FAILED")
        if not isinstance(go_2.defn, str):
            raise Exception("def TEST FAILED")
    else:
        if hasattr(go_1, 'defn'):
            raise Exception("def TEST FAILED")
        if hasattr(go_2, 'defn'):
            raise Exception("def TEST FAILED")
    out.write("GOTerm: {}\n".format(go_1.__repr__()))
    out.write("GOTerm: {}\n".format(go_2.__repr__()))

#################################################################
# main
#################################################################
if __name__ == '__main__':
    test_all()


