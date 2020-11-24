import { Classes, Position } from "@blueprintjs/core";
import React, { FC } from "react";
import { LeftAlignedDetailsCell } from "../common/style";
import { FieldValues, StyledButton, StyledPopover } from "./style";

interface Props {
  values: string[];
}

const Popover: FC<Props> = ({ values }) => {
  return (
    <LeftAlignedDetailsCell>
      <FieldValues>
        {values[0]}
        <br />
        {values[1]}
      </FieldValues>
      {values.length > 2 && (
        <StyledPopover
          boundary="window"
          modifiers={{
            hide: { enabled: false },
            preventOverflow: { enabled: false },
          }}
          position={Position.BOTTOM}
          popoverClassName={Classes.POPOVER_CONTENT_SIZING}
          content={
            <div>
              {values.map((val, idx) => (
                <React.Fragment key={val}>
                  {val}
                  {idx !== values.length - 1 && <br />}
                </React.Fragment>
              ))}
            </div>
          }
        >
          <StyledButton minimal>+{values.length - 2}</StyledButton>
        </StyledPopover>
      )}
    </LeftAlignedDetailsCell>
  );
};

export default Popover;
