import { Button, Menu, MenuItem, Popover, Position } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { EXTERNAL_LINKS } from "src/common/constants/routes";
import { MenuWrapper } from "./style";

const LearnButton = (): JSX.Element => {
  return (
    <Popover content={<Content />} position={Position.BOTTOM}>
      <Button
        minimal
        rightIcon={IconNames.CHEVRON_DOWN}
        text="Help & Documentation"
      />
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
