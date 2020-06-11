using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using CSV, DataFrames, Pipe

languagesF = "../data/languages.csv"
formsF = "../data/forms.csv"

!(isfile(languagesF) && isfile(formsF)) && begin
      asjpZip = download(
            "https://zenodo.org/record/3843469/files/lexibank/asjp-v19.1.zip",
            "tmp/asjp-19.1.zip",
      )
      run(`unzip -o $asjpZip -d tmp/`)
      cp("tmp/lexibank-asjp-0c18d44/cldf/forms.csv", formsF)
      cp("tmp/lexibank-asjp-0c18d44/cldf/languages.csv", languagesF)
      for f in readdir("tmp")
            rm("tmp/" * f, recursive = true)
      end
end

forms = CSV.read(formsF)
languages = CSV.read(languagesF)

function cleanASJP(word)
      @pipe word |>
            replace(_, r"[ \*~\"]" => "") |>
            replace(_, r"(.)(.)(.)\$" => s"\2")
end

languages[!, :longname] = @pipe languages |>
      zip(_.classification_wals, _.Name) |>
      join.(_, ".") |>
      replace.(_, "-" => "_")

forms[!, :simplified] = cleanASJP.(forms.Value)

languages = languages[.!occursin.("Oth.", languages.classification_wals), :]

conceptCoverage = @pipe forms |>
      unique(_, [:Language_ID, :gloss_in_source]) |>
      groupby(_, :gloss_in_source) |>
      combine(nrow, _) |>
      sort(_, :nrow, rev = true)

concepts = conceptCoverage.gloss_in_source[1:40]


forms40 = forms[map(x -> x ∈ concepts, forms.gloss_in_source), :]

languageCoverage = @pipe forms40 |>
      unique(_, [:Language_ID, :gloss_in_source]) |>
      groupby(_, :Language_ID) |>
      combine(nrow, _) |>
      sort(_, :nrow, rev = true)

doculects = languageCoverage.Language_ID[languageCoverage.nrow.>=30]

forms40 = forms40[map(x -> x ∈ doculects, forms40.Language_ID), :]

asjpLong = innerjoin(
      forms40,
      languages[:, [:Name, :longname]],
      on = :Language_ID => :Name,
)

asjpWide = @pipe asjpLong |>
      groupby(_, [:longname, :gloss_in_source]) |>
      combine(:simplified => x -> join(x, "-"), _) |>
      unstack(_, :longname, :gloss_in_source, :simplified_function)

CSV.write("../data/asjp19wide.csv", asjpWide)
