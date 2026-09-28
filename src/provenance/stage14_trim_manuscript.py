from pathlib import Path

root = Path(__file__).resolve().parent
body = root / 'paper/sections/body_reviewfriendly.tex'
s = body.read_text(encoding='utf-8')

start = s.index('\\section{Noise consequence}')
end = s.index('\\section{Discussion}', start)
s = s[:start] + s[end:]

intro_old = 'The deterministic switch is the object of the theorem and finite certificates; the noise experiment asks how such a switch appears in reported frequency estimates.'
s = s.replace(intro_old, 'The deterministic switch is the object of the theorem and finite certificates.')
s = s.replace('Classical frequency-estimation thresholds and ambiguity are established \\cite{rife1974,williamson1994,xu2004}. Accordingly, the noisy experiment below is a consequence of this deterministic boundary, not a new universal threshold law. ', '')
s = s.replace('Sec.~S6 summarizes the full replay.', 'Sec.~S5 summarizes the full replay.')
start = s.index('Near deterministic equal cost, noise chooses')
end = s.index('\\section{Limitations}', start)
s = s[:start] + 'The continuum theorem shows that the mechanism persists for every sufficiently large odd record, although its proof does not give an explicit threshold. The centered midpoint kernel fixes the first two large-record corrections to the crossing coefficient. It does not certify a connected branch-exchange boundary at finite normalized spacing.\n\n' + s[end:]
s = s.replace('The phase cusp, resolution sweep, asymmetries, and noisy fits have the numerical scope', 'The phase cusp, resolution sweep, and asymmetries have the numerical scope')
body.write_text(s, encoding='utf-8')

evidence = root / 'paper/sections/evidence_table.tex'
s = evidence.read_text(encoding='utf-8')
s = '\n'.join(line for line in s.splitlines() if not line.startswith(('AWGN branch selection', 'Frozen-projector branch-selection'))) + '\n'
evidence.write_text(s, encoding='utf-8')

supp = root / 'supplement/supplement.tex'
s = supp.read_text(encoding='utf-8').replace('\\input{S5_noise_search_details.tex}\n\\input{S6_reproducibility.tex}', '\\input{S5_independent_replay.tex}')
supp.write_text(s, encoding='utf-8')

replay = root / 'supplement/S5_independent_replay.tex'
source = (root / 'supplement/S6_reproducibility.tex').read_text(encoding='utf-8')
replay.write_text(source.replace('\\section*{S6.', '\\section*{S5.'), encoding='utf-8')
