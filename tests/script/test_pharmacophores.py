from six import StringIO

from kripodb.script.pharmacophores import get_run


def test_get__onefrag():
    out_file = StringIO()

    get_run('data/pharmacophores.h5', '3wsj_MK1_frag1', out_file)

    result = out_file.getvalue()
    nr_lines = len(result.split('\n'))
    assert '3wsj_MK1_frag1' in result and \
           'LIPO' in result and \
           '$$$$' in result and \
           nr_lines == 50


def test_get__all():
    out_file = StringIO()

    get_run('data/pharmacophores.h5', None, out_file)

    result = out_file.getvalue()
    nr_lines = len(result.split('\n'))
    assert '3wsj_MK1_frag1' in result and \
           'LIPO' in result and \
           '$$$$' in result and \
           nr_lines == 17611
