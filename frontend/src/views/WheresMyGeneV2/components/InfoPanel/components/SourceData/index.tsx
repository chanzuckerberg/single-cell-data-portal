import { List, ListItem } from "@czi-sds/components";
import { useContext, useMemo } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  aggregateCollectionsFromDatasets,
  useFilterDimensions,
} from "src/common/queries/wheresMyGene";
import { StateContext } from "src/views/WheresMyGeneV2/common/store";
import { Content, InfoText, ListSubheader, Wrapper } from "./style";
import { ROUTES } from "src/common/constants/routes";

interface Collection {
  name: string;
  url: string;
  datasets: { id: string; label: string }[];
}

interface Collections {
  [name: string]: Collection;
}

export default function SourceData(): JSX.Element {
  const { selectedFilters } = useContext(StateContext);
  const { data: filterDimensions } = useFilterDimensions();
  let { datasets = [] } = filterDimensions;
  if (selectedFilters.datasets.length > 0)
    datasets = datasets.filter((dataset) =>
      selectedFilters.datasets.includes(dataset.id)
    );
  const collections: Collections = useMemo(() => {
    return aggregateCollectionsFromDatasets(datasets);
  }, [datasets]);

  return (
    <Wrapper>
      <Content>
        <InfoText>
          Gene Expression is powered by primary data submitted to CZ CELLxGENE
          Discover. See exceptions and processing notes{" "}
          <a
            href={ROUTES.WMG_DOCS_DATA_PROCESSING}
            target="_blank"
            rel="noreferrer noopener"
          >
            here
          </a>
          .
        </InfoText>
        <List ordered data-testid="source-data-list">
          {Object.values(collections).map(({ name, url, datasets }) => {
            return (
              <ListItem ordered key={name} fontSize="xxxs">
                <List
                  key={name}
                  subheader={
                    <a
                      href={url}
                      target="_blank"
                      rel="noopener"
                      onClick={() => {
                        track(EVENTS.VIEW_COLLECTION_PAGE_CLICKED, {
                          collection_name: name,
                        });
                      }}
                    >
                      <ListSubheader>{name}</ListSubheader>
                    </a>
                  }
                >
                  {datasets.map(({ id, label }) => {
                    return (
                      <ListItem fontSize="xxxs" key={id}>
                        {label}
                      </ListItem>
                    );
                  })}
                </List>
              </ListItem>
            );
          })}
        </List>
      </Content>
    </Wrapper>
  );
}
