import {
  Button,
  Intent,
  Menu,
  MenuItem,
  Popover,
  Position,
} from "@blueprintjs/core";
import { EXTERNAL_LINKS } from "src/common/constants/routes";
import { MenuWrapper } from "./style";

const LearnButton = (): JSX.Element => {
  return (
    <Popover position={Position.BOTTOM} content={<Content />}>
      <Button minimal intent={Intent.PRIMARY}>
        Help & Documentation
      </Button>
    </Popover>
  );
};

function Content() {
  return (
    <MenuWrapper>
      <Menu>
        <MenuItem
          data-test-id="docs-data-portal"
          text="Documentation"
          href={EXTERNAL_LINKS.DOCS_DATA_PORTAL}
          target="_blank"
          rel="noopener"
        />
        <MenuItem
          data-test-id="docs-tutorials"
          text="Tutorials"
          href={EXTERNAL_LINKS.DOCS_TUTORIAL}
          target="_blank"
          rel="noopener"
        />
        <MenuItem
          data-test-id="docs-roadmap"
          text="Our Roadmap"
          href={EXTERNAL_LINKS.DOCS_ROADMAP}
          target="_blank"
          rel="noopener"
        />
      </Menu>
    </MenuWrapper>
  );
}

export default LearnButton;
