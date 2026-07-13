"""
Dengue Genome Aggregator - Single-Sequence Analysis Web App

Local Flask app: paste/upload one DENV sequence (complete genome or partial
fragment), get back a full computational-virology workup rendered as an
illustrated report. See sequence_analysis.py for the pipeline itself.

Ephemeral by design: nothing submitted here is written to disk or to
dengue_genomes.db (opened read-only throughout -- see sequence_analysis's
PRAGMA query_only connections). Local research tool only; no auth, binds to
127.0.0.1 by default.

Usage:
    python app.py
"""

import logging
import random
import sqlite3
import sys
from pathlib import Path

from flask import Flask, jsonify, render_template, request
from werkzeug.exceptions import RequestEntityTooLarge

import fingerprints
from sequence_analysis import DB_PATH, SequenceInputError, analyze_sequence

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# A single DENV genome is ~11KB, so this cap is not "how big is one genome" --
# sequence_analysis._clean_input only ever fully parses the first FASTA record
# regardless of file size, so this is purely a backstop against someone
# pointing the upload box at something enormous by mistake (e.g. one of the
# 100MB+ per-serotype Nextstrain exports in dengue_data/).
MAX_INPUT_BYTES = 50_000_000

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = MAX_INPUT_BYTES + 1_000  # + slack for JSON wrapper overhead


@app.errorhandler(RequestEntityTooLarge)
def handle_too_large(_e):
    return jsonify({'error': (
        f'Request body exceeds the {MAX_INPUT_BYTES:,}-byte limit. A single DENV genome is only '
        f'~11KB, so this is likely a multi-genome batch FASTA file -- trim it down to the one '
        f'record you want and re-upload.'
    )}), 413


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/api/analyze', methods=['POST'])
def api_analyze():
    payload = request.get_json(silent=True) or {}
    sequence = payload.get('sequence', '')

    if not isinstance(sequence, str) or not sequence.strip():
        return jsonify({'error': 'No sequence provided.'}), 400

    if len(sequence.encode('utf-8', errors='ignore')) > MAX_INPUT_BYTES:
        return jsonify({'error': (
            f'Input is larger than {MAX_INPUT_BYTES:,} bytes. A single DENV genome is only '
            f'~11KB, so this is likely a multi-genome batch FASTA file rather than one sequence -- '
            f'only the first record would be analyzed anyway, so trim the file down to the one '
            f'record you want and re-upload.'
        )}), 400

    try:
        result = analyze_sequence(sequence)
    except SequenceInputError as e:
        return jsonify({'error': str(e)}), 400
    except Exception:
        logger.exception("Unexpected error analyzing sequence")
        return jsonify({'error': 'Unexpected server error while analyzing this sequence.'}), 500

    return jsonify(result)


@app.route('/api/example')
def api_example():
    try:
        serotype = int(request.args.get('serotype', 0))
    except ValueError:
        serotype = 0
    if serotype not in (1, 2, 3, 4):
        return jsonify({'error': 'serotype must be 1, 2, 3, or 4.'}), 400

    conn = sqlite3.connect(str(DB_PATH))
    conn.execute("PRAGMA query_only = TRUE")
    try:
        rows = conn.execute(
            "SELECT accession, sequence, country, collection_date FROM genomes "
            "WHERE serotype = ? AND sequence IS NOT NULL AND length(sequence) > 10000",
            (serotype,),
        ).fetchall()
    finally:
        conn.close()

    if not rows:
        return jsonify({'error': f'No example genomes found for serotype {serotype}.'}), 404

    accession, sequence, country, collection_date = random.choice(rows)
    return jsonify({
        'accession': accession,
        'sequence': sequence,
        'country': country,
        'collection_date': collection_date,
        'fasta': f'>{accession} (example DENV-{serotype} genome, {country or "unknown origin"})\n{sequence}',
    })


def _check_startup_requirements():
    if not DB_PATH.exists():
        logger.error(f"Database not found at {DB_PATH}. Run main.py first to build it.")
        sys.exit(1)
    if not fingerprints.all_caches_exist():
        logger.error(
            "Fingerprint cache not found in fingerprint_cache/. "
            "Run `python build_fingerprints.py` first, then restart this app."
        )
        sys.exit(1)


if __name__ == '__main__':
    _check_startup_requirements()
    logger.info(f"Using database: {DB_PATH}")
    app.run(host='127.0.0.1', port=5000, debug=True)
