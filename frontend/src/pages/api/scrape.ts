import type { NextApiRequest, NextApiResponse } from "next";
import cheerio from "cheerio";

let cellTypeToWikipediaUrl: { [cellType: string]: string } = {};

class SPARQLQueryDispatcher {
  endpoint: string;
  constructor(endpoint: string) {
    this.endpoint = endpoint;
  }

  query(sparqlQuery: string) {
    const fullUrl = this.endpoint + "?query=" + encodeURIComponent(sparqlQuery);
    const headers = { Accept: "application/sparql-results+json" };

    return fetch(fullUrl, { headers }).then((body) => body.json());
  }
}

const endpointUrl = "https://query.wikidata.org/sparql";
const sparqlQuery = `SELECT ?item  ?sitelink1 ?sitelinks ?clid
WHERE 
{
    ?item wdt:P279* wd:Q7868 .  
    FILTER NOT EXISTS {?item wdt:P279* wd:Q729463 .} # Excluding sensory receptors, as Wikidata mixes cells and other non-uniquely-cell receptors.

{ ?sitelink1 schema:isPartOf <https://en.wikipedia.org/>;
    schema:about ?item;}

            {?item wdt:P7963 ?clid }

?item rdfs:label ?label
FILTER (lang(?label) = 'en')
        
?item wikibase:sitelinks ?sitelinks . 
        
}
ORDER BY DESC(?sitelinks)`;

const queryDispatcher = new SPARQLQueryDispatcher(endpointUrl);

async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    if (Object.keys(cellTypeToWikipediaUrl).length === 0) {
      const data = await queryDispatcher.query(sparqlQuery);
      const bindings = data.results.bindings;
      bindings.forEach((binding: any) => {
        const cellTypeId = binding.clid.value as string;
        const wikiUrl = binding.sitelink1.value as string;
        cellTypeToWikipediaUrl[cellTypeId] = wikiUrl;
      });
    }
    const cellTypeId = req.query.cellTypeId as string;
    const wikipediaUrl = cellTypeToWikipediaUrl[cellTypeId];
    if (!wikipediaUrl) {
      res.status(200).json({ content: "Try OLS" });
    } else {
      const response = await fetch(wikipediaUrl);
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }

      const html = await response.text();
      const $ = cheerio.load(html);

      const content = $(".mw-parser-output");
      const firstHeading = content.find("h2").first();
      const paragraphs = content
        .children()
        .slice(0, firstHeading.index() + 1)
        .filter("p");

      let text = paragraphs
        .map((_, el) => $(el).text().trim())
        .get()
        .join("\n\n");
      if (text.startsWith("\n\n")) text = text.slice(2);
      const cleanedText = text.replace(/\[\d+\]/g, "");
      res.status(200).json({
        content: cleanedText !== "" ? cleanedText : "Parsing error.",
      });
    }
  } catch (error) {
    console.error(`Error fetching the Wikipedia page: ${error}`);
    res.status(500).json({ error: "Error fetching the Wikipedia page" });
  }
}

export default handler;
