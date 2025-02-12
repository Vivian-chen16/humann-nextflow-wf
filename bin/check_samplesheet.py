#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FQ_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
    )

    VALID_TYPE_FORMATS = (
        "PE",
        "PEI",
        "SE"
    )

    def __init__(
        self,
        sample_col="sample",
        first_col="f1",
        second_col="f2",
        type_col="type",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the first (or only)
                FASTQ file path (default "fastq_1").
            second_col (str): The name of the column that contains the second (if any)
                FASTQ file path (default "fastq_2").
            type_col (str): The column that records whether the sample contains single-
                or paired-end sequencing reads or a genome (default "type").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._second_col = second_col
        self._type_col = type_col
        self._seen = set()
        self._seen_names = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_first(row)
        self._validate_second(row)
        self._validate_pair(row)
        self._validate_type(row)
        self._seen.add((row[self._sample_col], row[self._first_col]))
        self._seen_names.add((row[self._sample_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].strip().replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the first FASTQ/FASTA entry is non-empty and has the right format."""
        if len(row[self._first_col]) <= 0:
            raise AssertionError("At least the first FASTA/FASTQ file is required.")
        else:
            self._validate_fastq_format(row[self._first_col])

    def _validate_second(self, row):
        """Assert that the second FASTQ entry has the right format if it exists."""
        if len(row[self._second_col]) > 0:
            self._validate_fastq_format(row[self._second_col])

    def _validate_type(self, row):
        """Assert that the type entry is non-empty and has the right format."""
        if len(row[self._type_col]) <= 0:
            raise AssertionError("A value for the type field is required.")
        row[self._type_col] = row[self._type_col].upper()
        if row[self._type_col] not in self.VALID_TYPE_FORMATS:
            raise AssertionError(
                f"The type column requires a valid format. "
                f"Valid formats include (case-insensitive): "
                f"{', '.join(self.VALID_TYPE_FORMATS)}"
            )

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        if row[self._first_col] and row[self._second_col]:
            if row[self._type_col] != "PE":
                raise AssertionError("The type must be 'PE' if f1 and f2 are provided.")
            else:
                first_col_suffix = Path(row[self._first_col]).suffixes[-2:]
                second_col_suffix = Path(row[self._second_col]).suffixes[-2:]
                if first_col_suffix != second_col_suffix:
                    raise AssertionError("FASTQ pairs must have the same file extensions.")
        elif row[self._type_col] == "PE":
            raise AssertionError("If the type is 'PE' both f1 and f2 must be provided.")

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_FQ_FORMATS):
            raise AssertionError(
                f"The FASTQ file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_FQ_FORMATS)}"
            )
        
    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTA/FASTQ filename is unique.
        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of sample name and FASTA/FASTQ must be unique.")

    def validate_unique_samplenames(self):
        """
        Assert that sample names are all unique.
        """
        if len(self._seen_names) != len(self.modified):
            raise AssertionError("Sample names must be unique.")


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure:

            sample,f1,f2,type,background
            SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz,PE
            SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz,PE
            SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,SE
            SAMPLE_PEI,SAMPE_PEI_RUN1_1.fastq.gz,PEI
    """
    required_columns = {"sample", "f1", "f2", "type"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
        checker.validate_unique_samplenames()
    header = list(reader.fieldnames)

    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
