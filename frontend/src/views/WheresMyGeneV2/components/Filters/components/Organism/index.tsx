import { EMPTY_ARRAY } from "src/common/constants/utils";
import { StyledDropdown, Wrapper, Label } from "../common/style";
import { InputDropdownProps, Props, tempOnChange } from "./types";
import { useConnect } from "./connect";
import Dialog from "./components/Dialog";

export default function Organism({ isLoading }: Props): JSX.Element {
  const {
    handleDialogCancel,
    handleDialogConfirm,
    handleOnChange,
    isDialogOpen,
    organism,
    organisms,
    DropdownMenuProps,
  } = useConnect();

  return (
    <Wrapper>
      <Label>Organism</Label>
      <StyledDropdown
        label={organism?.name || "Select"}
        options={organisms || EMPTY_ARRAY}
        onChange={handleOnChange as tempOnChange}
        InputDropdownProps={{ ...InputDropdownProps, disabled: isLoading }}
        data-testid="add-organism"
        DropdownMenuProps={DropdownMenuProps}
      />
      <Dialog
        isOpen={isDialogOpen}
        handleCancel={handleDialogCancel}
        handleConfirm={handleDialogConfirm}
      />
    </Wrapper>
  );
}
