from graphviz import Digraph
import paramiko
import sys

class RemoteReader:
    def __init__(self):
        self.transport = paramiko.Transport('192.168.1.23', 22)
        self.transport.connect()
        self.transport.auth_password('jifahu', '3RtT32qc')
        self.sftp = self.transport.open_sftp_client()

    def open(self, remote):
        try:
            file = self.sftp.open(remote, 'r')
            return file
        except Exception as e:
            print(e)


if __name__ == '__main__':
    g = Digraph('Tree', filename='btree.gz', format='png')
    rr = RemoteReader()
    with rr.open('/vol1/agis/ruanjue_group/jifahu/project/pgrc_learn/build/tree.txt') as f:
        for line in f:
            s = line.strip().split()
            print(s)
            pointer, n = int(s[1]), int(s[2])
            if s[0] == 'external':
                keys = s[3:3+n]
                g.node(str(pointer), '%d:%s' % (pointer, ' '.join(keys)), color='lightgrey', stype='filled', shape='rectangle')
            else:
                keys, children = s[3:3+n], s[3+n:]
                g.node(str(pointer), '%d:%s' % (pointer, ' '.join(keys)), shape='rectangle')
                for c in children:
                    g.edge(str(pointer), str(c))
        g.render()
