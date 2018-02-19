from six import StringIO

from kripodb.script.pharmacophores import get_run, sd2phar


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


def test_sd2phar(example1_sdfile, example1_pharblock):
    out_file = StringIO()
    frag_id = 'some_frag_id'
    sd2phar(example1_sdfile, out_file, frag_id)

    assert out_file.getvalue() == example1_pharblock
