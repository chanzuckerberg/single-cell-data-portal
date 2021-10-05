import { Button, Intent, Menu, MenuItem, Popover } from "@blueprintjs/core";
import { EXTERNAL_LINKS } from "src/common/constants/routes";

const LearnButton = (): JSX.Element => {
  return (
    <Popover content={<Content />}>
      <Button minimal intent={Intent.PRIMARY}>
        Learn
      </Button>
    </Popover>
  );
};

function Content() {
  return (
    <Menu>
      <MenuItem
        data-test-id="docs-data-portal"
        text="Documentation"
        href={EXTERNAL_LINKS.DOCS_DATA_PORTAL}
        target="_blank"
        rel="noopener"
      />
      <MenuItem
        data-test-id="docs-roadmap"
        text="Roadmap"
        href={EXTERNAL_LINKS.DOCS_ROADMAP}
        target="_blank"
        rel="noopener"
      />
    </Menu>
  );
}

export default LearnButton;
