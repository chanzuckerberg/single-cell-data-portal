import { Button, Menu, MenuItem, Popover } from "@blueprintjs/core";

interface Props {
  onSelect: (genesetIndex: number) => void;
}

export default function Genesets({ onSelect }: Props): JSX.Element {
  function handleSelect(genesetIndex: number): void {
    onSelect(genesetIndex);
  }

  return (
    <Popover
      content={
        <Menu>
          <MenuItem onClick={() => handleSelect(0)} text="1" />
          <MenuItem onClick={() => handleSelect(1)} text="2" />
          <MenuItem onClick={() => handleSelect(2)} text="3" />
        </Menu>
      }
    >
      <Button>Gene Sets</Button>
    </Popover>
  );
}
