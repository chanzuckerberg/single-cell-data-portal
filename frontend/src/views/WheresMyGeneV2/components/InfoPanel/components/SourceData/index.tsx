import { List, ListItem } from "@czi-sds/components";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { Content, InfoText, ListSubheader, Wrapper } from "./style";
import { ROUTES } from "src/common/constants/routes";
import { useConnect } from "./connect";

export default function SourceData(): JSX.Element {
  const { collections } = useConnect();

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
