import os,sys


def main():
    for filename in sys.argv[1:]:
        f = open( filename )
        tokens = filename.split('.')
        tokens[-1] = 'modified.fa'
        out = open( '.'.join(tokens),'w')
        for r in f:
            if r.startswith('>'):
                r = r.replace(':','_')
                r = r.replace('-','_')
            out.write(r)
        f.close()
        out.close()



if __name__=='__main__':
    main()
