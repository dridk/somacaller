import argparse
import sys
from somacaller.model import SomaModel

class VafyParser(object):
    def __init__(self):

        parser = argparse.ArgumentParser(
            description="Vafy tools",
            usage="""vafy <command> [<args>] 
subcommand are : 
    fit               learn model from a bam list
    test           compute prediction score from a bam file 
    plot              plot model and prediction

            """,
        )

        parser.add_argument("command", help="subcommand to run")

        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command + "_cmd"):
            print("Unrecognized command")
            parser.print_help()
            exit(1)

        getattr(self, args.command + "_cmd")(sys.argv[2:])

    def fit_cmd(self, argv):
        parser = argparse.ArgumentParser(description="learn models")
        parser.add_argument('files', nargs='+')
        parser.add_argument('-r', '--reference', help="ex: hg19.fa", required=True)
        parser.add_argument('-t', '--target', help="ex: hotspot.bed", required=True)
        parser.add_argument('-o', '--output', help="ex: model.h5", required=True)

        args = parser.parse_args(argv)       

        print(args.files)

        model = SomaModel(args.target, args.reference)
        model.fit(args.files)
        model.to_hdf(args.output)   

    def test_cmd(self, argv):
        parser = argparse.ArgumentParser(description="test bam")
        parser.add_argument('file')
        parser.add_argument('-m', '--model', help="ex: model.h5", required=True)
        parser.add_argument('-o', '--output', help="ex: result.txt", required=True)
        parser.add_argument('-t', '--type', default="tab", choices=["tab","json"], help="output type")
        args = parser.parse_args(argv)       


        model = SomaModel.from_hdf(args.model)
        df = model.test(args.file)

        if args.type == "tab":
            df.to_csv(args.output, index=None)
        if args.type == "json":
            df.to_json(args.output,orient="records")

    def plot_cmd(self, argv):
        parser = argparse.ArgumentParser(description="test bam")
        parser.add_argument('-p', '--position', help="ex: chr3:23424", required=True)
        parser.add_argument('-m', '--model', help="ex: model.h5", required=True)
        parser.add_argument('-o', '--output', help="ex: graph.json", required=True)
        args = parser.parse_args(argv)   

        model = SomaModel.from_hdf(args.model)
        chart = model.plot(args.position)

        chart.save(args.output)



if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.DEBUG)
    VafyParser()