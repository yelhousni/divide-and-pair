package curve25519

import (
	"crypto/rand"
	"testing"
)

func TestFilippoAgreement(t *testing.T) {
	params := curveParameters()
	for i := 0; i < 20; i++ {
		k, _ := rand.Int(rand.Reader, &params.Order)
		var p PointAffine
		p.ScalarMultiplication(&params.Base, k)

		pornin := p.isInSubGroupPornin()
		porninF := p.isInSubGroupPorninFilippo()
		quarticExpF := p.isInSubGroupQuarticExpFilippo()
		quarticF := p.isInSubGroupQuarticFilippo()

		if !pornin || !porninF || !quarticExpF || !quarticF {
			t.Fatalf("subgroup: pornin=%v porninFilippo=%v quarticExpFilippo=%v quarticFilippo=%v", pornin, porninF, quarticExpF, quarticF)
		}
	}
	t.Log("all filippo variants agree")
}
